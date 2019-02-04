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

def create_window():
    """creates application window and contains all window functionality"""

    def restart():
        """Restarts window"""
        root.destroy()
        create_window()
        return

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






    # Run window indefinitely
    root.mainloop()



def main():
    # Main function to run.
    create_window()


if __name__ == '__main__':
    main()
