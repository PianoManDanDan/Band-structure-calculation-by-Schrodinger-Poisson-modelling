"""
This is a module for calculating the strain between consecutive layers
of grown semiconductor material. If the strain is greater than 5%,
an alert box pops up warning the user that the chosen materials used
may be difficult, or even impossible, to grow.
"""

__version__ = '0.1'
__author__ = 'Daniel Martin'
__all__ = ['strain']

import numpy as np

try:
    from tkinter import messagebox
except ImportError:
    from TKinter import tkMessageBox as messagebox


def strain(material_list):
    """
    Calculates strain between 2 layers of semiconductor material.

    Parameters:
    -----------
    material_list: list of size N
                   List of material classes in layers
    """

    for i in range(1, len(material_list)):
        mat_strain = material_list[i].lattice_constant / \
                     material_list[i-1].lattice_constant
        if mat_strain > 1.05 or mat_strain < 0.95:
            messagebox.showwarning("Warning - Lattice mismatch found",
                                   "Lattice mismatch greater than 5% "
                                   "found between {0} and {1}. This "
                                   "material configuration may not be "
                                   "possible to grow".format(
                                                material_list[i-1].name,
                                                material_list[i].name))
            break

    return


if __name__ == '__main__':
    import materials

    mat_list = [materials.AlGaAs(0.2, 10), materials.InSb(10),
                materials.InP(10), materials.AlGaAs(0.2, 10)]

    strain(mat_list)
