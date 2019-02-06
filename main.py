"""
This script makes use of the schrodinger_poisson module and 
is the main script to be run for the Schrodinger-Poisson solver. 

This script, when run, will create a GUI interface for the Schrodinger-
Poisson solver and run the necessary calculations.
"""
__author__ = 'Daniel Martin'
__version__ = '0.1'

# Imports
try:
    import tkinter as tk
except ImportError:
    import Tkinter as tk

import numpy as np
import matplotlib.pyplot as plt

# WORK OUT HOW TO GET IMPORTS FROM __INIT__.PY
from schrodinger_poisson.schrodinger import solve_schrodinger
# from schrodinger_poisson.poisson import solve_poisson
# from schrodinger_poisson.diffusion import solve_diffusion
# from schrodinger_poisson import schrodinger_poisson as SP


def create_window():
    """creates application window and contains all window functionality"""

    def restart():
        """Restarts window"""
        root.destroy()
        create_window()
        return

    def go():
        """Functionality for GO button in tkinter window"""
        # JSON material variables
        # x, V = file.open(path, unpack=True)
        N = N_states.get()
        print(N)

    # Window setup
    root = tk.Tk()
    root.title("Band Structure Calculator")
    # root.state('zoomed')
    root.wm_iconbitmap('icon.ico')

    # Add top menu toolbar
    menu = tk.Menu(root)
    root.config(menu=menu)
    file_menu = tk.Menu(menu, tearoff=False)
    menu.add_cascade(label="File", menu=file_menu)
    file_menu.add_command(label='New calculation...', command=restart)
    file_menu.add_separator()
    file_menu.add_command(label='Exit', command=root.destroy)

    # Left side options

    # N states to find
    N_states = tk.Spinbox(root, from_=1, to=10)  # Is 10 enough?
    N_states.pack()

    go_button = tk.Button(root, command=go, text='GO!')
    go_button.pack()

    # embedded matplotlib window


    # Run window indefinitely
    root.mainloop()


def main():
    # Main function to run.
    create_window()


if __name__ == '__main__':
    main()
