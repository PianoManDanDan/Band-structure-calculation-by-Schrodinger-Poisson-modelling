"""
This script makes use of the schrodinger_poisson module and 
is the main script to be run for the Schrodinger-Poisson solver. 

This script, when run, will create a GUI interface for the Schrodinger-
Poisson solver and run the necessary calculations.
"""
__author__ = 'Daniel Martin'
__version__ = '0.1'

# Imports
import warnings

try:
    import tkinter as tk
    from tkinter import messagebox
    from tkinter import ttk
except ImportError:
    import Tkinter as tk
    from TKinter import tkMessageBox as messagebox
    from Tkinter import ttk

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk
from scipy import constants

import schrodinger_poisson as sp


class Window:
    """class object for creating tkinter window"""

    def __init__(self):
        # VARIABLE SETUP
        # These variables are assigned values in the go() function.
        # Initialised in __init__ for better readability.
        self.material_list = None
        self.thickness = None
        self.total_thickness = None
        self.me = None
        self.mhh = None
        self.mlh = None
        self.Eg = None
        self.dielectric = None
        self.x = None
        self.V = None
        self.T = None
        self.eigvals = None
        self.eigvects = None
        self.N_points = None
        self.N_states = None
        self.particle = None
        self.num_materials = 1

        self.reset_settings()

        # TKINTER SETUP
        self.root = tk.Tk()
        self.root.grid_rowconfigure(0, weight=1)
        self.root.grid_columnconfigure(1, weight=1)

        self.root.title('Band Structure Calculator')
        # self.root.state('zoomed') # Uncomment to open fullscreen
        self.root.wm_iconbitmap('icon.ico')

        # Add top menu toolbar
        self.menu = tk.Menu(self.root)
        self.root.config(menu=self.menu)
        self.file_menu = tk.Menu(self.menu, tearoff=False)
        self.menu.add_cascade(label='File', menu=self.file_menu)
        self.file_menu.add_command(label='New calculation...',
                                   command=self.restart)
        self.file_menu.add_separator()
        self.file_menu.add_command(label='Exit',
                                   command=self.root.destroy)

        self.options_menu = tk.Menu(self.menu, tearoff=False)
        self.menu.add_cascade(label='Settings', menu=self.options_menu)
        self.options_menu.add_cascade(label='Adjust Settings',
                                      command=self.adjust_settings)
        self.options_menu.add_cascade(label='Reset Defaults',
                                      command=self.reset_settings)

        # LEFT HAND SIDE
        self.LHS = tk.Frame(self.root)
        self.LHS.grid(row=0, column=0, sticky='news')
        self.vcmd = (self.LHS.register(self.validate_entry), '%P')
        # Temperature box
        temp_frame = tk.Frame(self.LHS)
        temp_frame.grid(row=0, column=0, columnspan=7)
        temperature_label = tk.Label(temp_frame, text='Temperature (K)')
        temperature_label.grid(row=0, column=0, padx=10, pady=5,
                               sticky='n')
        self.temperature_entry = tk.Entry(temp_frame, validate='key',
                                          validatecommand=self.vcmd)
        self.temperature_entry.bind('<KeyRelease>',
                                    self.set_temperature)
        self.temperature_entry.grid(row=0, column=1)

        # Labels for adding in materials
        material_label = tk.Label(self.LHS, text='Material')
        material_label.grid(row=1, column=0, pady=10)
        x_label = tk.Label(self.LHS, text='x')
        x_label.grid(row=1, column=1, pady=10)
        thickness_label = tk.Label(self.LHS, text='Thickness\n(nm)')
        thickness_label.grid(row=1, column=2, pady=10)
        Eg_label = tk.Label(self.LHS, text='Eg\n(eV)')
        Eg_label.grid(row=1, column=3, pady=10)
        me_label = tk.Label(self.LHS, text='me\n(m0)')
        me_label.grid(row=1, column=4, pady=10)
        mh_label = tk.Label(self.LHS, text='mh\n(m0)')
        mh_label.grid(row=1, column=5, pady=10)
        mlh_label = tk.Label(self.LHS, text='mlh\n(m0)')
        mlh_label.grid(row=1, column=6, pady=10)
        dielectric_label = tk.Label(self.LHS, text=u'\u03f5')
        dielectric_label.grid(row=1, column=7, pady=10)

        # Set material parameters
        self.material_choices = ['InSb', 'GaAs',
                                 u'Al\u2093Ga\u208d\u2081\u208b\u2093\u208eAs',
                                 'InP', 'InAs', 'GaSb', 'Custom']

        self.materials_dropdown = [ttk.Combobox(self.LHS,
                                                state='disabled',
                                                values=self.material_choices,
                                                width=10)]
        self.materials_dropdown[0].bind('<<ComboboxSelected>>',
                                        self.insert_parameters)
        self.x_entry = [tk.Entry(self.LHS, validate='key',
                                 state='disabled',
                                 validatecommand=self.vcmd, width=10)]
        self.x_entry[0].bind('<KeyRelease>', self.insert_parameters)
        self.thickness_entry = [tk.Entry(self.LHS, state='disabled',
                                         validate='key',
                                         validatecommand=self.vcmd,
                                         width=10)]
        self.Eg_entry = [tk.Entry(self.LHS, state='disabled',
                                  validate='key',
                                  validatecommand=self.vcmd, width=10)]
        self.me_entry = [tk.Entry(self.LHS, state='disabled',
                                  validate='key',
                                  validatecommand=self.vcmd, width=10)]
        self.mh_entry = [tk.Entry(self.LHS, state='disabled',
                                  validate='key',
                                  validatecommand=self.vcmd, width=10)]
        self.mlh_entry = [tk.Entry(self.LHS, state='disabled',
                                   validate='key',
                                   validatecommand=self.vcmd, width=10)]
        self.dielectric_entry = [tk.Entry(self.LHS, state='disabled',
                                          validate='key',
                                          validatecommand=self.vcmd,
                                          width=10)]

        self.materials_dropdown[0].grid(row=2, column=0,
                                        padx=(5, 2), pady=(4, 0))
        self.x_entry[0].grid(row=2, column=1,
                             padx=2, pady=(4, 0))
        self.thickness_entry[0].grid(row=2, column=2,
                                     padx=2, pady=(4, 0))
        self.Eg_entry[0].grid(row=2, column=3,
                              padx=2, pady=(4, 0))
        self.me_entry[0].grid(row=2, column=4,
                              padx=2, pady=(4, 0))
        self.mh_entry[0].grid(row=2, column=5,
                              padx=2, pady=(4, 0))
        self.mlh_entry[0].grid(row=2, column=6,
                               padx=2, pady=(4, 0))
        self.dielectric_entry[0].grid(row=2, column=7,
                                      padx=(2, 5), pady=(4, 0))

        # Button to add in new layers
        self.add_layer_button = tk.Button(self.LHS,
                                          command=self.add_layer,
                                          text='+', padx=5)
        self.add_layer_button.grid(row=3, column=0, pady=5)
        self.remove_layer_button = tk.Button(self.LHS,
                                             command=self.remove_layer,
                                             text='-',
                                             state='disabled', padx=8)
        self.remove_layer_button.grid(row=3, column=7, pady=5)

        # Calculate Potential Button and GO! button
        self.calc_pot_button = tk.Button(self.LHS,
                                         command=self.calculate_potential,
                                         text='Calculate Potential',
                                         state='disabled',
                                         padx=5, pady=5)
        self.calc_pot_button.grid(row=4, column=0, columnspan=4)
        self.go_button = tk.Button(self.LHS, command=self.go,
                                   text='GO!', state='disabled',
                                   padx=5, pady=5)
        self.go_button.grid(row=4, column=4, columnspan=4)

        self.RHS = tk.Frame(self.root)
        self.RHS.grid(row=0, column=1, sticky='news')

        self.RHS.grid_rowconfigure(1, weight=1)
        self.RHS.grid_columnconfigure(0, weight=1)

        # matplotlib window
        fig = plt.Figure(figsize=(5.3, 4.3))
        self.figure = fig.add_subplot(111)
        self.figure.set_xlabel('Growth axis (nm)')
        self.figure.set_ylabel('E (eV)')
        self.figure.grid()
        self.plot = FigureCanvasTkAgg(fig, master=self.RHS)
        self.plot.get_tk_widget().grid(row=1, column=0, sticky='nesw',
                                       padx=10, pady=(0, 10))
        self.toolbar_frame = tk.Frame(master=self.RHS)
        self.toolbar_frame.grid(row=0, column=0, padx=10, pady=(10, 0),
                                sticky='news')
        self.toolbar = NavigationToolbar2Tk(self.plot,
                                            self.toolbar_frame)
        self.toolbar.update()
        self.toolbar.grid(row=0, column=0, sticky='news')
        self.plot.draw()

        # mainloop
        self.root.mainloop()

    def restart(self):
        """Restarts window"""
        self.root.destroy()
        self.__init__()
        return

    def reset_settings(self):
        """Reset N_points and N_states to default settings"""
        self.N_points = 1000
        self.N_states = 5
        self.particle = 'e'
        return

    def adjust_settings(self):
        """
        Adjusts user setting for the number of points and the number of
        states to find.
        """

        def accept():
            """Functionality for accept button"""
            try:
                self.N_states = int(N_states_spinbox.get())
                self.particle = particle.get()
                self.N_points = int(N_points_box.get())
            except ValueError:
                self.N_points = 1000
                warnings.warn('Number of points not specified. '
                              'Defaulting to 1000 points',
                              SyntaxWarning, stacklevel=1)
            print(particle.get())
            options_window.destroy()
            return

        def validate_int(value_if_allowed):
            """Validates entry into N_points_box is an integer"""
            if value_if_allowed == '':
                return True

            try:
                int(value_if_allowed)
                return True
            except ValueError:
                return False

        options_window = tk.Tk()
        options_window.title('Options')
        options_window.wm_iconbitmap('icon.ico')
        vcmd = (options_window.register(validate_int), '%P')
        tk.Label(options_window, text='No. of points:',
                 borderwidth=5).grid(row=0, column=0, pady=5)
        N_points_box = tk.Entry(options_window, validate='key',
                                validatecommand=vcmd)
        N_points_box.grid(row=0, column=1, padx=5, sticky='ew')

        tk.Label(options_window, text='No. of states to find:',
                 borderwidth=5).grid(row=1, column=0, pady=5)
        N_states_default = tk.StringVar(options_window)
        N_states_default.set(str(self.N_states))
        N_states_spinbox = tk.Spinbox(options_window, from_=1,
                                      to=self.N_points,
                                      textvariable=N_states_default,
                                      state='readonly')
        N_states_spinbox.grid(row=1, column=1, padx=5, sticky='ew')

        particle = tk.StringVar(options_window)
        particle.set('e')
        tk.Label(options_window, text='Particle:',
                 borderwidth=5).grid(row=2, column=0, rowspan=3,
                                     pady=5, sticky='n')
        e_radio = tk.Radiobutton(options_window, text='Electron',
                                 variable=particle, value='e')
        e_radio.grid(row=2, column=1, pady=(5, 0), sticky='w')
        mhh_radio = tk.Radiobutton(options_window, text='Heavy Hole',
                                   variable=particle, value='hh')
        mhh_radio.grid(row=3, column=1, sticky='w')
        mlh_radio = tk.Radiobutton(options_window, text='Light Hole',
                                   variable=particle, value='lh')
        mlh_radio.grid(row=4, column=1, pady=(0, 5), sticky='w')

        accept_button = tk.Button(options_window, command=accept,
                                  text='Accept', padx=5, pady=5)
        accept_button.grid(row=5, column=0, columnspan=2,
                           sticky='ns', padx=5, pady=8)

        return

    @staticmethod
    def validate_entry(value_if_allowed):
        """
        Allows only numerical characters to be entered into entry boxes
        """
        if value_if_allowed == '':
            return True

        try:
            float(value_if_allowed)
            return True
        except ValueError:
            return False

    def set_temperature(self, event):
        """sets self.T as the user types in a value"""
        if self.temperature_entry.get() == '' or \
                float(self.temperature_entry.get()) <= 0:
            self.materials_dropdown[-1].configure(state='disabled')
            self.x_entry[-1].configure(state='readonly')
            self.thickness_entry[-1].configure(state='readonly')
            self.Eg_entry[-1].configure(state='readonly')
            self.me_entry[-1].configure(state='readonly')
            self.mh_entry[-1].configure(state='readonly')
            self.mlh_entry[-1].configure(state='readonly')
            self.dielectric_entry[-1].configure(state='readonly')
            return
        else:
            self.T = float(self.temperature_entry.get())
            self.materials_dropdown[-1].configure(state='readonly')
            if self.materials_dropdown[-1].get() == \
                    u'Al\u2093Ga\u208d\u2081\u208b\u2093\u208eAs':
                self.x_entry[-1].configure(state='normal')
            self.thickness_entry[-1].configure(state='normal')
            self.Eg_entry[-1].configure(state='normal')
            self.me_entry[-1].configure(state='normal')
            self.mh_entry[-1].configure(state='normal')
            self.mlh_entry[-1].configure(state='normal')
            self.dielectric_entry[-1].configure(state='normal')

        # Automatically update Entry boxes if self.T changed
        for i in range(self.num_materials):
            material = self.materials_dropdown[i].get()
            if material == '' or material == 'Custom' or \
                    material == u'Al\u2093Ga\u208d\u2081\u208b\u2093\u208eAs':
                continue
            else:
                self.Eg_entry[i].configure(state='normal')
                self.Eg_entry[i].delete(0, 'end')
                material = sp.materials[material](self.T)
                self.Eg_entry[i].insert(0, str(material.Eg)[:5])
                self.Eg_entry[i].configure(state='readonly')
            self.Eg_entry[-1].configure(state='normal')
        return

    def add_layer(self):
        """
        Adds functionality to '+' button. Adds a new material layer
        row to the window. Currently allows up to a maximum of 10
        material layers.
        """
        # Only add if previous layer is filled
        if self.materials_dropdown[-1].get() == '' \
                or self.thickness_entry[-1].get() == '':
            return

        # Remove '+' button from window
        self.add_layer_button.grid_remove()
        self.remove_layer_button.grid_remove()
        self.calc_pot_button.grid_remove()
        self.go_button.grid_remove()

        # disable previous layer
        self.materials_dropdown[-1].configure(state='disabled')
        self.x_entry[-1].configure(state='readonly')
        self.thickness_entry[-1].configure(state='readonly')
        self.Eg_entry[-1].configure(state='readonly')
        self.me_entry[-1].configure(state='readonly')
        self.mh_entry[-1].configure(state='readonly')
        self.mlh_entry[-1].configure(state='readonly')
        self.dielectric_entry[-1].configure(state='readonly')

        # Add layer
        self.materials_dropdown.append(ttk.Combobox(self.LHS,
                                                    state='readonly',
                                                    values=self.material_choices,
                                                    width=10))
        self.materials_dropdown[-1].bind('<<ComboboxSelected>>',
                                         self.insert_parameters)
        self.x_entry.append(tk.Entry(self.LHS, validate='key',
                                     state='readonly',
                                     validatecommand=self.vcmd,
                                     width=10))
        self.x_entry[-1].bind('<KeyRelease>', self.insert_parameters)
        self.thickness_entry.append(tk.Entry(self.LHS, validate='key',
                                             validatecommand=self.vcmd,
                                             width=10))
        self.Eg_entry.append(tk.Entry(self.LHS, validate='key',
                                      validatecommand=self.vcmd,
                                      width=10))
        self.me_entry.append(tk.Entry(self.LHS, validate='key',
                                      validatecommand=self.vcmd,
                                      width=10))
        self.mh_entry.append(tk.Entry(self.LHS, validate='key',
                                      validatecommand=self.vcmd,
                                      width=10))
        self.mlh_entry.append(tk.Entry(self.LHS, validate='key',
                                       validatecommand=self.vcmd,
                                       width=10))
        self.dielectric_entry.append(tk.Entry(self.LHS, validate='key',
                                              validatecommand=self.vcmd,
                                              width=10))

        if self.num_materials < 25:
            self.materials_dropdown[-1].grid(row=self.num_materials + 2,
                                             column=0, padx=(5, 2),
                                             pady=(4, 0))
            self.x_entry[-1].grid(row=self.num_materials + 2, column=1,
                                  padx=2, pady=(4, 0))
            self.thickness_entry[-1].grid(row=self.num_materials + 2,
                                          column=2, padx=2, pady=(4, 0))
            self.Eg_entry[-1].grid(row=self.num_materials + 2, column=3,
                                   padx=2, pady=(4, 0))
            self.me_entry[-1].grid(row=self.num_materials + 2, column=4,
                                   padx=2, pady=(4, 0))
            self.mh_entry[-1].grid(row=self.num_materials + 2, column=5,
                                   padx=2, pady=(4, 0))
            self.mlh_entry[-1].grid(row=self.num_materials + 2,
                                    column=6, padx=2, pady=(4, 0))
            self.dielectric_entry[-1].grid(row=self.num_materials + 2,
                                           column=7, padx=(2, 5),
                                           pady=(4, 0))

            self.remove_layer_button.grid(row=self.num_materials + 3,
                                          column=7, pady=(10, 0))
            self.remove_layer_button.configure(state='normal')
            self.calc_pot_button.grid(row=self.num_materials + 4,
                                      column=0, columnspan=4,
                                      padx=5, pady=(10, 0))
            self.go_button.grid(row=self.num_materials + 4, column=4,
                                columnspan=4, padx=5, pady=(10, 0))
        if self.num_materials < 24:
            self.add_layer_button.grid(row=self.num_materials + 3,
                                       column=0, pady=(10, 0))

        self.num_materials += 1

        if self.num_materials > 2:
            self.calc_pot_button.configure(state='normal')
        return

    def remove_layer(self):
        """
        Adds functionality to '-' button. Removes last material layer
        row in the window.
        """
        # Don't remove last layer
        if self.num_materials < 2:
            return

        # Remove buttons from window
        self.add_layer_button.grid_remove()
        self.remove_layer_button.grid_remove()
        self.calc_pot_button.grid_remove()
        self.go_button.grid_remove()

        # Remove last layer
        self.materials_dropdown[-1].grid_remove()
        self.x_entry[-1].grid_remove()
        self.thickness_entry[-1].grid_remove()
        self.Eg_entry[-1].grid_remove()
        self.me_entry[-1].grid_remove()
        self.mh_entry[-1].grid_remove()
        self.mlh_entry[-1].grid_remove()
        self.dielectric_entry[-1].grid_remove()

        self.materials_dropdown.pop(-1)
        self.x_entry.pop(-1)
        self.thickness_entry.pop(-1)
        self.Eg_entry.pop(-1)
        self.me_entry.pop(-1)
        self.mh_entry.pop(-1)
        self.mlh_entry.pop(-1)
        self.dielectric_entry.pop(-1)

        # Reactivate new final layer
        self.materials_dropdown[-1].configure(state='readonly')
        if self.materials_dropdown[-1].get() == \
                u'Al\u2093Ga\u208d\u2081\u208b\u2093\u208eAs':
            self.x_entry[-1].configure(state='normal')
        self.thickness_entry[-1].configure(state='normal')
        self.Eg_entry[-1].configure(state='normal')
        self.me_entry[-1].configure(state='normal')
        self.mh_entry[-1].configure(state='normal')
        self.mlh_entry[-1].configure(state='normal')
        self.dielectric_entry[-1].configure(state='normal')

        # re-add buttons
        self.add_layer_button.grid(row=self.num_materials+2,
                                   column=0, pady=(10, 0))
        self.calc_pot_button.grid(row=self.num_materials + 3, column=0,
                                  columnspan=4, padx=5, pady=(10, 0))
        self.go_button.grid(row=self.num_materials + 3, column=4,
                            columnspan=4, padx=5, pady=(10, 0))
        if self.num_materials > 1:
            self.remove_layer_button.grid(row=self.num_materials+2,
                                          column=7, pady=(10, 0))


        self.num_materials -= 1

        if self.num_materials < 2:
            self.remove_layer_button.configure(state='disabled')

        if self.num_materials < 3:
            self.calc_pot_button.configure(state='disabled')
            self.go_button.configure(state='disabled')
        return

    def insert_parameters(self, event):
        """
        Inserts parameters into entry boxes when material is selected.
        """

        chosen_material = self.materials_dropdown[-1].get()

        # Clear values in case of material change
        self.Eg_entry[-1].delete(0, 'end')
        self.me_entry[-1].delete(0, 'end')
        self.mh_entry[-1].delete(0, 'end')
        self.mlh_entry[-1].delete(0, 'end')
        self.dielectric_entry[-1].delete(0, 'end')

        if chosen_material == '' or chosen_material == 'Custom':
            return
        if chosen_material == \
                u'Al\u2093Ga\u208d\u2081\u208b\u2093\u208eAs' and \
                self.x_entry[-1].get() == '':
            self.x_entry[-1].configure(state='normal')
            return
        if chosen_material == \
                u'Al\u2093Ga\u208d\u2081\u208b\u2093\u208eAs':
            x_entry = float(self.x_entry[-1].get())
            material = sp.materials['AlGaAs'](x_entry, self.T)
        else:
            self.x_entry[-1].delete(0, 'end')
            self.x_entry[-1].configure(state='disabled')
            material = sp.materials[chosen_material](self.T)

        self.Eg_entry[-1].insert(0, str(material.Eg)[:5])
        self.me_entry[-1].insert(0, str(material.me)[:5])
        self.mh_entry[-1].insert(0, str(material.mhh)[:5])
        self.mlh_entry[-1].insert(0, str(material.mlh)[:5])
        self.dielectric_entry[-1].insert(0, str(material.dielectric_constant)[:5])

        return

    def set_parameters(self):
        """set parameters into variables for use in calculations"""

        # Make list of material classes
        self.material_list = []
        for i in range(self.num_materials):
            material = self.materials_dropdown[i].get()
            if material == \
                    u'Al\u2093Ga\u208d\u2081\u208b\u2093\u208eAs':
                x_entry = float(self.x_entry[i].get())
                self.material_list.append(sp.materials['AlGaAs'](
                                                            x_entry,
                                                            self.T))
            elif material == 'Custom':
                self.material_list.append(sp.materials[material](
                    float(self.Eg_entry[i].get()),
                    float(self.me_entry[i].get()),
                    float(self.mh_entry[i].get()),
                    float(self.mlh_entry[i].get()),
                    float(self.dielectric_entry[i].get()))
                )
            else:
                self.material_list.append(sp.materials[material](self.T))

        # Make array of thicknesses
        self.thickness = np.array([])
        for i in range(self.num_materials):
            self.thickness = np.append(self.thickness,
                                       float(self.thickness_entry[i].get()))
        self.thickness *= 1e-9
        self.total_thickness = sum(self.thickness)
        cum_thickness = np.insert(self.thickness, 0, 0)
        cum_thickness = np.cumsum(cum_thickness)

        # create growth axis parameter
        self.x = np.linspace(0, self.total_thickness, self.N_points)

        # Make array of Eg
        self.Eg = np.zeros_like(self.x)
        for i in range(self.num_materials):
            self.Eg[self.x >= cum_thickness[i]] = float(self.Eg_entry[i].get())
        self.Eg *= constants.eV

        # Make array of electron mass
        self.me = np.zeros_like(self.x)
        for i in range(self.num_materials):
            self.me[self.x >= cum_thickness[i]] = float(self.me_entry[i].get())
        self.me *= constants.m_e

        # Make array of heavy hole mass
        self.mhh = np.zeros_like(self.x)
        for i in range(self.num_materials):
            self.mhh[self.x >= cum_thickness[i]] = float(self.mh_entry[i].get())
        self.mhh *= constants.m_e

        # Make array of light hole mass
        self.mlh = np.zeros_like(self.x)
        for i in range(self.num_materials):
            self.mlh[self.x >= cum_thickness[i]] = float(self.mlh_entry[i].get())
        self.mlh *= constants.m_e

        # Make array of dielectric
        self.dielectric = np.zeros_like(self.x)
        for i in range(self.num_materials):
            self.mhh[self.x >= cum_thickness[i]] = float(self.dielectric_entry[i].get())
        self.dielectric *= constants.epsilon_0

        return

    def calculate_potential(self):
        """Calculates potential from input materials"""

        # Remove last level if row is empty
        if self.materials_dropdown[-1].get() == '':
            self.remove_layer()

        # Check all values filled in if last layer is 'Custom'
        if self.materials_dropdown[-1].get() == 'Custom':
            params = np.array([self.Eg_entry[-1].get(),
                               self.me_entry[-1].get(),
                               self.mh_entry[-1].get(),
                               self.mlh_entry[-1].get(),
                               self.dielectric_entry[-1].get()])
            if any(params == ''):
                raise ValueError('Final layer must have a full set of '
                                 'parameters')

        # Check for final layer thickness
        if self.thickness_entry[-1].get() == '':
            raise ValueError('Final layer must have a thickness')

        # Set parameters from entry boxes
        self.set_parameters()

        # Strain calculation
        sp.strain(self.material_list)

        # calculate V
        self.V = sp.anderson(self.x, self.material_list,
                             self.thickness)

        # Plot potential to plotting window.
        self.figure.cla()
        self.figure.plot(self.x / 1e-9, self.V / constants.eV, 'k')
        self.figure.set_xlabel('x (nm)')
        self.figure.set_ylabel('E (eV)')
        self.figure.grid()
        self.plot.draw()


        # Enable Go button
        self.go_button.configure(state='normal')
        return

    def go(self):
        """Functionality for GO! button. Calculates band structure"""

        # Determine what particle to use: electron, heavy hole,
        # light hole
        particle_options = {'e': self.me, 'hh': self.mhh, 'lh': self.mlh}
        particle = particle_options[self.particle]

        # Calculate eigenvalues and eigenvectors for potential
        self.eigvals, self.eigvects = sp.solve_schrodinger(self.x,
                                                           self.V,
                                                           particle)

        # Plot eigenvectors
        num_plotted = 0
        i = 0
        while num_plotted < self.N_states:
            if self.eigvals[i] > max(self.V):
                break
            if self.eigvals[i] < 0:
                i += 1
                continue
            self.figure.plot(self.x / 1e-9,
                             self.eigvects[:, i] + self.eigvals[i] / constants.eV,
                             label=r'$\psi_{0}$'.format(num_plotted))
            num_plotted += 1
            i += 1
        self.figure.legend(loc='best')
        self.plot.draw()

        return


if __name__ == '__main__':
    Window()
