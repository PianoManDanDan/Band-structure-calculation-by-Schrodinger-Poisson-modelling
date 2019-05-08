"""
This is a module for calculating the potential wells between layers
of grown semiconductor material using Anderson's Rule.
"""

__version__ = '0.1'
__author__ = 'Daniel Martin'
__all__ = ['anderson', 'valence_band']

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

    V = np.zeros_like(growth_axis)

    for i in range(len(dV)):
        V[growth_axis >= cum_thickness[i]] = material_list[0].Eg + \
                                              dV[i]


    V *= constants.eV
    return V


def valence_band(conduction_band, growth_axis, material_list,
                 material_thickness):
    """ - Redo!!!
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
    conduction_band = np.copy(conduction_band)
    conduction_band /= constants.eV

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
            Eg[i] = Eg[i-1] * 0.66
        elif prev_mat.name == 'GaAs' and next_mat.name == 'AlGaAs':
            Eg[i] = Eg[i-1] / 0.66
        else:
            pass

    dV = np.asarray([Eg[i] - Eg[0] for i in range(len(Eg))])

    V2 = np.zeros_like(conduction_band)
    for i in range(len(Eg)):
        V2[growth_axis >= cum_thickness[i]] = \
            conduction_band[growth_axis >= cum_thickness[i]][0] - \
            2 * Eg[i]
        V2[growth_axis >= cum_thickness[i]] = -material_list[0].Eg - \
                                             dV[i]

    V2 *= constants.eV
    return V2


if __name__ == '__main__':
    import materials
    import matplotlib.pyplot as plt
    from schrodinger import solve_schrodinger

    mat_list = [materials.AlGaAs(0.2, 10), materials.GaAs(10),
                materials.InP(10), materials.AlGaAs(0.2, 10)]
    # mat_list = [materials.GaAs(10), materials.InSb(10),
    #             materials.InP(10), materials.GaAs(10)]
    # mat_list = [materials.AlGaAs(0.2, 10), materials.GaAs(10),
    #               materials.GaAs(10), materials.AlGaAs(0.2, 10)]
    thickness = np.array([10, 10, 10, 10]) * 1e-9

    x = np.linspace(0, sum(thickness), 1000)
    V = anderson(x, mat_list, thickness)
    V2 = valence_band(V, x, mat_list, thickness)

    m = np.ones_like(x) * (0.063 + 0.083 * 0.2)
    m[x >= 10e-9] = 0.063
    m[x >= 10e-9] = 0.08
    m[x >= 30e-9] = (0.063 + 0.083 * 0.2)
    m *= constants.m_e
    e1, e2 = solve_schrodinger(x, -V2, m)

    print(e1[:5] / constants.eV)

    plt.figure()
    plt.plot(x / 1e-9, V / constants.eV, 'k')
    plt.plot(x / 1e-9, V2 / constants.eV, 'k')
    plt.plot(x/1e-9, -e1[:5]/constants.eV - e2[:,:5])
    plt.grid()
    plt.show()

