"""
This is a module for calculating the potential wells between layers
of grown semiconductor material using Anderson's Rule.
"""

__version__ = '0.1'
__author__ = 'Daniel Martin'
__all__ = ['anderson']

import numpy as np
from scipy import constants


def anderson(growth_axis, material_list, material_thickness):
    """
    Calculates potential based on Anderson's Rule. If adjacent layers
    are GaAs/AlGaAs, gives back 66/33% relationship between band gaps.

    Parameters:
    -----------
    growth_axis: array of size N
                 linearly spaced points along growth axis.
    material_list: list of size M
                   List of material classes in layer
    material_thickness: array of size M
                        Array of thickness values for each material
                        layer (m).

    Out:
    ----
    V: Array of size N
       Potential profile from materials due to Anderson's Rule
    """

    material_thickness = np.insert(material_thickness, 0, 0)
    cum_thickness = np.cumsum(material_thickness)
    dV = np.ones(len(material_thickness)-1) * material_list[0].Eg / 2

    for i in range(1, len(material_list)):
        prev_mat = material_list[i-1]
        next_mat = material_list[i]
        if prev_mat.name == 'AlGaAs' and next_mat.name == 'GaAs':
            dV[i] = 0.34 * dV[i-1]
        elif prev_mat.name == 'GaAs' and next_mat.name == 'AlGaAs':
            dV[i] = dV[i-1] / 0.34
        else:
            dV[i] = dV[i-1] + (next_mat.Eg - prev_mat.Eg) / 2

    V = np.zeros_like(growth_axis)
    for i in range(len(dV)):
        V[growth_axis >= cum_thickness[i]] = dV[i]

    if material_list[0].name == material_list[-1].name:
        V[growth_axis >= cum_thickness[-2]] = V[0]

    V *= constants.eV
    return V


if __name__ == '__main__':
    import materials
    import matplotlib.pyplot as plt

    mat_list = [materials.AlGaAs(0.2, 10), materials.GaAs(10),
                materials.InP(10), materials.AlGaAs(0.2, 10)]
    # mat_list = [materials.GaAs(10), materials.InSb(10),
    #             materials.InP(10), materials.GaAs(10)]
    # mat_list = [materials.AlGaAs(0.2, 10), materials.GaAs(10),
    #               materials.GaAs(10), materials.AlGaAs(0.2, 10)]
    thickness = np.array([10, 10, 10, 10]) * 1e-9

    x = np.linspace(0, sum(thickness), 1000)
    V = anderson(x, mat_list, thickness)

    plt.figure()
    plt.plot(x, V, 'k')
    plt.grid()
    plt.show()

