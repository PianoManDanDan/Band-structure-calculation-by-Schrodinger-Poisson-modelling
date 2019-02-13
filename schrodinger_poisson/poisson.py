"""
This is a module for solving the Poisson equation. This module is not
currently implemented into the main program.
"""

__version__ = '0.1'
__author__ = 'Daniel Martin'
__all__ = ['charge_density', 'E_field', 'potential_deformation',
           'solve_poisson']

import numpy as np
from scipy import constants


def charge_density(q, x, wavefunction, volume_density):
    """
    Calculates the spatial charge density for a static charge being
    added into the potential profile.
    """
    if len(x) < 2:
        raise ValueError("Spatial coordinates must be more than a "
                         "single point")
    if len(x) != len(volume_density):
        raise ValueError("Spatial coordinates and volume density must "
                         "have same shape")

    N = len(x)
    dz = abs(x[1] - x[0])

    return q * (N * wavefunction**2 - volume_density) * dz


def E_field(charge_density_, material_permittivity):
    """
    Calculates the electric field produced by placing a charge within
    a potential well.
    """
    return charge_density_ / (2 * material_permittivity)


def potential_deformation():
    """
    Calculates the deformation of the potential due to the addition of
    the charge.
    """
    pass


def solve_poisson(x, V, q=1.6e-19):
    """
    Combines all other functions in module to output a new potential
    profile.
    """
    return None


if __name__ == '__main__':
    import matplotlib.pyplot as plt

    x = np.linspace(-1e-8, 1e-8, 1000)
    V = np.ones_like(x) * 200e-3 * 1.6e-19
    V[abs(x) < 5e-9] = 0

    V_new = solve_poisson(x, V)

    plt.figure()
    plt.plot(x, V_new)
    plt.show()




