"""
This script makes use of the schrodinger_poisson module and 
is the main script to be run for the Schrodinger-Poisson solver. 

This script, when run, will create a GUI interface for the Schrodinger-
Poisson solver and run the necessary calculations.
"""
__author__ = 'Daniel Martin'
__version__ = '0.1'

# Imports
import json
import warnings

try:
    import tkinter as tk
    from tkinter.filedialog import askopenfilename
    from tkinter import ttk
except ImportError:
    import Tkinter as tk
    from Tkinter.tkFileDialog import askopenfilename
    from Tkinter import ttk

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk

import schrodinger_poisson as sp


class Window:
    """class object for creating tkinter window"""

    def __init__(self):
        # VARIABLE SETUP
        # These variables are assigned values in the go() function.
        # Initialised in __init__ for better readability.
        self.me = None
        self.mhh = None
        self.mlh = None
        self.Eg = None
        self.dielectric_constant = None
        self.x = None
        self.V = None
        self.T = None
        self.N_points = None
        self.N_states = None
        self.num_materials = 3

        self.reset_settings()

        # TKINTER SETUP
        self.root = tk.Tk()
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


        #####################################################################################
        # Set material parameters
        self.material_choices = ['InSb', 'GaAs', 'AlGaAs', 'InP',
                                 'InAs', 'GaSb', 'Custom']
        self.material_dropdown1 = ttk.Combobox(self.LHS,
                                               state='readonly',
                                               values=self.material_choices,
                                               width=10)
        self.material_dropdown1.grid(row=2, column=0, padx=(5, 2))
        self.x_entry1 = tk.Entry(self.LHS, validate='key',
                                 state='disabled',
                                 validatecommand=self.vcmd, width=10)
        self.x_entry1.grid(row=2, column=1, padx=2, pady=2)
        self.thickness_entry = tk.Entry(self.LHS, validate='key',
                                        validatecommand=self.vcmd,
                                        width=10)
        self.thickness_entry.grid(row=2, column=2, padx=2, pady=2)
        self.Eg_entry = tk.Entry(self.LHS, validate='key',
                                 validatecommand=self.vcmd, width=10)
        self.Eg_entry.grid(row=2, column=3, padx=2, pady=2)
        self.me_entry = tk.Entry(self.LHS, validate='key',
                                 validatecommand=self.vcmd, width=10)
        self.me_entry.grid(row=2, column=4, padx=2, pady=2)
        self.mh_entry = tk.Entry(self.LHS, validate='key',
                                 validatecommand=self.vcmd, width=10)
        self.mh_entry.grid(row=2, column=5, padx=2, pady=2)
        self.mlh_entry = tk.Entry(self.LHS, validate='key',
                                  validatecommand=self.vcmd, width=10)
        self.mlh_entry.grid(row=2, column=6, padx=2, pady=2)
        self.dielectric_entry = tk.Entry(self.LHS, validate='key',
                                         validatecommand=self.vcmd,
                                         width=10)
        self.dielectric_entry.grid(row=2, column=7, padx=2, pady=2)
        #####################################################################################
        # Button to add in new layers
        self.add_layer_button = tk.Button(self.LHS,
                                          command=self.add_layer,
                                          text='+', padx=5)
        self.add_layer_button.grid(row=5, column=0, pady=5)

        """
        # choose material
        self.material = tk.StringVar(self.root)
        self.material.set('---')
        tk.Label(self.root, text='Material profile:', borderwidth=10) \
            .grid(row=0, column=0)
        material_choices = ['InSb', 'GaAs', 'AlGaAs', 'InP', 'InAs',
                            'GaSb', 'Select own material...']
        # self.material_dropdown = tk.OptionMenu(self.root, self.material,
        #                                        *material_choices,
        #                                        command=print('working'))
        self.material_dropdown = ttk.Combobox(self.root,
                                              state='readonly',
                                              values=material_choices)
        self.material_dropdown.grid(row=0, column=1, padx=10, pady=5)

        # potential select
        self.potential_button = tk.Button(self.root,
                                          command=self.get_potential,
                                          text='Calculate Potential',
                                          padx=10, pady=5)
        self.potential_button.grid(row=1, column=0, columnspan=2,
                                   padx=10, pady=5)

        # add go button
        self.go_button = tk.Button(self.root, command=self.go,
                                   text='GO!', padx=20, pady=10)
        self.go_button.grid(row=3, column=0, columnspan=2,
                            padx=10, pady=5)
        """
        self.RHS = tk.Frame(self.root)
        self.RHS.grid(row=0, column=1, sticky='news')

        # matplotlib window
        fig = plt.Figure(figsize=(5.3, 4.3))
        self.figure = fig.add_subplot(111)
        self.figure.format_coord = lambda x, y: ''
        self.figure.set_xlabel('Growth axis (nm)')
        self.figure.set_ylabel('E (meV)')
        self.figure.grid()
        self.plot = FigureCanvasTkAgg(fig, master=self.RHS)
        self.plot.get_tk_widget().grid(row=1, column=0, sticky='nesw',
                                       padx=10, pady=(0, 10))
        self.toolbar_frame = tk.Frame(master=self.RHS)
        self.toolbar_frame.grid(row=0, column=0, padx=10, pady=(10,0))
        self.toolbar = NavigationToolbar2Tk(self.plot, self.toolbar_frame)
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
        return

    def adjust_settings(self):

        def accept():
            try:
                self.N_points = int(N_points_box.get())
                self.N_states = int(N_states_spinbox.get())
            except ValueError:
                self.N_points = 1000
                warnings.warn('Number of points not specified. '
                              'Defaulting to 1000 points',
                              SyntaxWarning, stacklevel=1)
            options_window.destroy()
            return

        options_window = tk.Tk()
        options_window.title('Options')
        options_window.wm_iconbitmap('icon.ico')

        tk.Label(options_window, text='No. of points:',
                 borderwidth=5).grid(row=0, column=0, pady=5)
        N_points_box = tk.Entry(options_window)
        N_points_box.grid(row=0, column=1, padx=5)

        tk.Label(options_window, text='No. of states to find:',
                 borderwidth=5).grid(row=1, column=0, pady=5)
        N_states_spinbox = tk.Spinbox(options_window, from_=1,
                                      to=self.N_points)
        N_states_spinbox.grid(row=1, column=1, padx=5)

        accept_button =tk.Button(options_window, command=accept,
                                 text='Accept', padx=5, pady=5)
        accept_button.grid(row=2, column=0, columnspan=2,
                           sticky='ns', padx=5, pady=8)

        return

    @staticmethod
    def validate_entry(value_if_allowed):
        if value_if_allowed == '':
            return True

        try:
            float(value_if_allowed)
            return True
        except ValueError:
            return False

    def set_temperature(self, event):
        if self.temperature_entry.get() == '':
            return
        else:
            self.T = float(self.temperature_entry.get())
        return

    def add_layer(self):
        self.add_layer_button.grid_remove()
        label = tk.Label(self.LHS, text=str(self.num_materials))
        label.grid(row=self.num_materials+1, column=0)
        print(self.num_materials)
        if self.num_materials < 10:
            self.add_layer_button.grid(row=self.num_materials+2,
                                       column=0, pady=10)
        self.num_materials += 1
        return

    """
    def material_select(self):
        # Functionality for GO button in tkinter window
        # JSON material variables
        material = self.material.get()
        if material == 'Select own material...':
            material_path = askopenfilename(initialdir='./',
                                            title='Select Material File',
                                            filetypes=(
                                                ('JSON file', '*.json'),
                                                ('all files', '*.*')))
            with open(material_path) as f:
                material_data = json.load(f)
            self.me = material_data['me']
            self.mhh = material_data['mhh']
            self.mlh = material_data['mlh']
            self.Eg = material_data['Eg']
            self.dielectric_constant = \
                material_data['dielectric_constant']
        else:
            # material_x.get()  # material x box
            # material_T.get()  # material T box
            material_x = 0
            material_T = 0
            chosen_material = sp.materials.materials[material](
                material_x, material_T)
            self.me = chosen_material.me
            self.mhh = chosen_material.mhh
            self.mlh = chosen_material.mlh
            self.Eg = chosen_material.Eg
            self.dielectric_constant = \
                chosen_material.dielectric_constant
        return

    def get_potential(self):
        potential_path = askopenfilename(initialdir='./',
                                         title='Select Potential File',
                                         filetypes=(
                                             ('CSV file', '*.csv'),
                                             ('all files', '*.*')))
        x, V = np.loadtxt(potential_path, delimiter=',', unpack=True)
        self.x = x
        self.V = V

        # Plot x and V for sanity check - do in window!!
        # Update plot
        self.figure.cla()
        self.figure.plot(self.x/1e-9, self.V/1.6e-22, 'k')
        self.figure.set_xlabel('x (nm)')
        self.figure.set_ylabel('E (meV)')
        self.figure.grid()
        self.plot.draw()
        return

    def go(self):
        # Functionality for GO button in tkinter window
        if self.material.get() == '---':
            tk.messagebox.showwarning('Warning', 'Please select a '
                                                 'material')
            return
        elif np.any(self.V == None):
            tk.messagebox.showwarning('Warning', 'Please select a '
                                                 'potential profile')
            return
        self.material_select()

        N = self.N_states.get()
        print(N, self.Eg)

        return
    """

if __name__ == '__main__':
    Window()
