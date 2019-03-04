"""
This is a module for solving Schrodinger's equation using the matrix
method.

The main function of the module (solve_schrodinger) returns the energy
eigenvalues and corresponding wavefunctions for a potential, V, in a
domain, x.
"""

__version__ = '0.1'
__author__ = 'Daniel Martin'
__all__ = ['anderson']

import numpy as np


def anderson(growth_axis, material_list):
    """
    Calculates potential based on Anderson's Rule. If adjacent layers
    are GaAs/AlGaAs, gives back 66/33% relationship between band gaps.

    Parameters:
    -----------
    growth_axis: array of size N
                 linearly spaced points along growth axis.
    material_list: list
                   List of material classes in layer

    Out:
    ----
    V: Array of size N
       Potential profile from materials due to Anderson's Rule
    """

    # return V

if __name__ == '__main__':
    print('anderson')

