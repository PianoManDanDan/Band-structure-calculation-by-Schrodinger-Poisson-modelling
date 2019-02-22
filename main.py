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

# import
try:
    import tkinter as tk
    from tkinter.filedialog import askopenfilename
except ImportError:
    import Tkinter as tk
    from Tkinter.tkFileDialog import askopenfilename

import numpy as np
import matplotlib.pyplot as plt

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

        # TKINTER SETUP
        self.root = tk.Tk()
        self.root.title("Band Structure Calculator")
        # self.root.state('zoomed')
        self.root.wm_iconbitmap('icon.ico')

        # Add top menu toolbar
        self.menu = tk.Menu(self.root)
        self.root.config(menu=self.menu)
        self.file_menu = tk.Menu(self.menu, tearoff=False)
        self.menu.add_cascade(label="File", menu=self.file_menu)
        self.file_menu.add_command(label='New calculation...',
                                   command=self.restart)
        self.file_menu.add_separator()
        self.file_menu.add_command(label='Exit',
                                   command=self.root.destroy)

        # Left side options

        # choose material
        self.material = tk.StringVar(self.root)
        self.material.set('---')
        tk.Label(self.root, text='Material profile:', borderwidth=10) \
            .grid(row=0, column=0)
        material_choices = ['InSb', 'GaAs', 'AlGaAs', 'Si',
                            'Select own material...']
        # Add more choices with materials.py + Get Select own material
        # working
        self.material_dropdown = tk.OptionMenu(self.root, self.material,
                                               *material_choices)
        self.material_dropdown.grid(row=0, column=1, padx=10, pady=5)

        # potential select
        tk.Label(self.root, text='Potential profile:', borderwidth=10) \
            .grid(row=1, column=0)
        # self.potential_profile = ???
        # self.potential_profile.grid(row=1, column=1, padx=10, pady=5)

        # N states to find
        tk.Label(self.root, text='No. of states to find:',
                 borderwidth=10).grid(row=2, column=0)
        self.N_states = tk.Spinbox(self.root, from_=1, to=1000)
        self.N_states.grid(row=2, column=1, padx=10, pady=5)

        # add go button
        self.go_button = tk.Button(self.root, command=self.go,
                                   text='GO!', padx=20, pady=10)
        self.go_button.grid(row=3, column=0, columnspan=2,
                            padx=10, pady=5)

        # matplotlib window

        # mainloop
        self.root.mainloop()

    def restart(self):
        """Restarts window"""
        self.root.destroy()
        self.__init__()
        return

    def material_select(self):
        """Functionality for GO button in tkinter window"""
        # JSON material variables
        material = self.material.get()
        if material == '---':
            tk.messagebox.showwarning('Warning', 'Please select a '
                                                 'material')
            return
        elif material == 'Select own material...':
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
        return

    def go(self):
        """Functionality for GO button in tkinter window"""
        self.material_select()
        self.get_potential()

        N = self.N_states.get()
        print(N, self.Eg)

        return


def main():
    # Main function to run.
    Window()


if __name__ == '__main__':
    main()
