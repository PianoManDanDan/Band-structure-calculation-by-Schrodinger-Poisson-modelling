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

# from schrodinger_poisson.schrodinger import solve_schrodinger
# from schrodinger_poisson.poisson import solve_poisson
# from schrodinger_poisson.diffusion import solve_diffusion
# from schrodinger_poisson import schrodinger_poisson
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
        self.N_points = None
        self.N_states = None

        self.reset_settings()

        # TKINTER SETUP
        self.root = tk.Tk()
        self.root.title('Band Structure Calculator')
        # self.root.state('zoomed')
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
        # Left side options

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

        # matplotlib window
        fig = plt.Figure(figsize=(5.3, 4.3))
        self.figure = fig.add_subplot(111)
        self.figure.set_xlabel('x (nm)')
        self.figure.set_ylabel('E (meV)')
        self.figure.grid()
        self.plot = FigureCanvasTkAgg(fig, master=self.root)
        self.plot.get_tk_widget().grid(row=0, column=2, rowspan=10,
                                       columnspan=2, sticky='nesw',
                                       padx=10, pady=10)
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

    def material_select(self):
        """Functionality for GO button in tkinter window"""
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
        """Functionality for GO button in tkinter window"""
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


if __name__ == '__main__':
    Window()
