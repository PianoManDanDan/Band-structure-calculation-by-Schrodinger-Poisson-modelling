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
    """creates application window"""

    # Window setup
    root = tk.Tk()
    root.title("Band Structure Calculator")
    root.state('zoomed')
    root.wm_iconbitmap('icon.ico')



    # Run window indefinitely
    root.mainloop()



def main():
    # Main function to run.
    create_window()


if __name__ == '__main__':
    main()
