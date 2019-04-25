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

    ???
    V2: Array of size N
        Transition energy V...?
    """

    material_thickness = np.insert(material_thickness, 0, 0)
    cum_thickness = np.cumsum(material_thickness)
    Eg = [material_list[i].Eg for i in range(len(material_list))]
    Eg = np.asarray(Eg)

    for i in range(1, len(material_list)):
        prev_mat = material_list[i-1]
        next_mat = material_list[i]

        # Check if same material used twice in a row
        if next_mat.name == prev_mat.name:
            # Check if AlGaAs and check for different x values
            if next_mat.name == 'AlGaAs' and next_mat.x != prev_mat.x:
                continue

            Eg[i] = Eg[i-1]
            continue

        # 66%/33% split if GaAs-AlGaAs boundary
        if prev_mat.name == 'AlGaAs' and next_mat.name == 'GaAs':
            Eg[i] = Eg[i-1] * 0.33
        elif prev_mat.name == 'GaAs' and next_mat.name == 'AlGaAs':
            Eg[i] = Eg[i-1] / 0.33
        else:
            pass

    dV = np.asarray([Eg[i] - Eg[0] for i in range(len(Eg))])
    max_range = np.ptp(dV)

    # dV2 = max_range / ()

    V1 = np.zeros_like(growth_axis)
    # V2 = np.zeros_like(growth_axis)
    for i in range(len(dV)):
        V1[growth_axis >= cum_thickness[i]] = material_list[0].Eg + \
                                              dV[i]
        # V2[growth_axis >= cum_thickness[i]] = material_list[0].Eg + \
        #                                       dV[i] - \
        #                                       2 * material_list[i].Eg

    V1 *= constants.eV
    # V2 *= constants.eV
    return V1#, V2
    """UNCOMMENT V2 CODE IF TRANSITION ENERGIES BEING DONE!!!"""




if __name__ == '__main__':
    import materials
    import matplotlib.pyplot as plt
    import time

    mat_list = [materials.AlGaAs(0.2, 10), materials.GaAs(10),
                materials.InP(10), materials.AlGaAs(0.2, 10)]
    # mat_list = [materials.GaAs(10), materials.InSb(10),
    #             materials.InP(10), materials.GaAs(10)]
    # mat_list = [materials.AlGaAs(0.2, 10), materials.GaAs(10),
    #               materials.GaAs(10), materials.AlGaAs(0.2, 10)]
    thickness = np.array([10, 10, 10, 10]) * 1e-9

    x = np.linspace(0, sum(thickness), 1000)
    V = anderson(x, mat_list, thickness)
    # V1, V2 = anderson(x, mat_list, thickness)

    plt.figure()
    plt.plot(x / 1e-9, V / constants.eV, 'k')
    # plt.plot(x / 1e-9, V2 / constants.eV, 'k')
    plt.grid()
    plt.show()

