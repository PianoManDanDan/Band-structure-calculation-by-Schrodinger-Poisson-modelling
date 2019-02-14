"""
This is a module for solving the Poisson equation. This module is not
currently implemented into the main program.
"""

__version__ = '0.1'
__author__ = 'Daniel Martin'
__all__ = ['charge_density', 'E_field', 'material_permittivity',
           'potential_deformation', 'solve_poisson']

import numpy as np
from scipy import constants


def material_permittivity(relative_permittivity):
    """
    Calculates the overall permittivity of the material from the
    material's relative permittivity at a specific temperature.

    Parameters:
    -----------
    relative_permittivity: float
                           static dielectric constant for material at
                           chosen temperature.

    Out:
    ----
    permittivity: float
                  permittivity of material (F/m)
    """

    return relative_permittivity * constants.epsilon_0


def charge_density(x, wavefunction, volume_density, q):
    """
    Calculates the spatial charge density for a static charge being
    added into the potential profile.

    Parameters:
    -----------
    x: array of size N
       linearly spaced points defining spatial domain of well (m)
    wavefunction: array of size N
                  single wavefunction eigenstate of potential (units??)
    volume_density: array of size N
                    volume density of dopants at each position in
                    potential well (units??)
    q: float
           Charge of carriers being added to the system (C)

    Out:
    ----
    charge_density_: array of size N
                     net charge density of the semiconductor layer (
                     units??)  <-- is this an accurate description?
    """

    if len(x) < 2:
        raise ValueError('Spatial coordinates must be more than a '
                         'single point')
    if len(x) != len(volume_density):
        raise ValueError('Spatial coordinates and volume density must '
                         'have same shape')

    N = len(x)
    dz = abs(x[1] - x[0])

    return q * (N * wavefunction**2 - volume_density) * dz


def E_field(charge_density_, material_permittivity_):
    """
    Calculates the electric field produced by placing a charge within
    a potential well.

    Parameters:
    -----------
    charge_density_: array of size N
                     net charge density of the semiconductor layer (
                     units??)  <-- is this an accurate description?
    material_permittivity_: float
                            overall permittivity of material at given
                            temperature. If you have the relative
                            permittivity of the material, you can
                            call the material_permittivity() function
                            to calculate the overall permittivity.

    Out:
    ----
    E_field: array of size N
             Electric field along axis of charge density. (units??)
    """

    return charge_density_ / (2 * material_permittivity_)


def potential_deformation(x, E_field_):
    """
    Calculates the deformation of the potential due to the addition of
    the charge.

    Parameters:
    -----------
    x: array of size N
       linearly spaced points defining spatial domain of well (m)
    E_field: array of size N
             Electric field along axis of charge density. (units??)

    Out:
    ----
    poisson_potential: array of size N
                       additional potential caused by the inclusion of
                       a doping charge in the potential well.
    """

    return - np.trapz(E_field_, x)


def solve_poisson(x, V, wavefunction, volume_density,
                  relative_permittivity, q=1.6e-19):
    """
    Combines all other functions in module to output a new potential
    profile.

    Parameters:
    -----------
    x: array of size N
       linearly spaced points defining spatial domain of well (m)
    V: array of size N
       values of potential for each corresponding x value (J)
    wavefunction: array of size N
                  single wavefunction eigenstate of potential (units??)
    volume_density: array of size N
                    volume density of dopants at each position in
                    potential well (units??)
    relative_permittivity: float
                           static dielectric constant for material at
                           chosen temperature.
    q: float
       Charge of carriers being added to the system (C)

    Out:
    ----
    V_new: array of size N
           New potential to be put through Schrodinger solver
    """

    # DO VARIABLE TESTS

    permittivity = material_permittivity(relative_permittivity)
    charge_density_ = charge_density(q, x, wavefunction, volume_density)
    E_field_ = E_field(charge_density_, permittivity)
    deformation = potential_deformation(x, E_field_)
    return V + deformation


if __name__ == '__main__':
    import matplotlib.pyplot as plt

    x = np.linspace(-1e-8, 1e-8, 1000)
    V = np.ones_like(x) * 200e-3 * 1.6e-19
    V[abs(x) < 5e-9] = 0

    V_new = solve_poisson(x, V)

    plt.figure()
    plt.plot(x, V_new)
    plt.show()




