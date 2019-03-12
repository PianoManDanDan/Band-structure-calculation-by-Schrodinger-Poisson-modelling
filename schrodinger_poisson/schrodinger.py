"""
This is a module for solving Schrodinger's equation using the matrix
method. 

The main function of the module (solve_schrodinger) returns the energy
eigenvalues and corresponding wavefunctions for a potential, V, in a
domain, x.
"""

__version__ = '0.1'
__author__ = 'Daniel Martin'
__all__ = ['_beta', '_M', 'energy_eigenvalues', 'wavefunctions',
           'solve_schrodinger']

import numpy as np
from numpy import linalg as LA
from scipy import constants


def _beta(x, V, m, hbar):
    """
    This function calculates the beta value matrix.

    Parameters:
    -----------
    x: array of size N
       linearly spaced points defining spatial domain of well (m)
    V: array of size N
       values of potential for each corresponding x value (J)
    m: array of size N
       electron mass (kg)
    hbar: float
          Planck's constant / 2 pi (Js)

    Out:
    ----
    beta: array of size N
          array of beta values to go into calculating M
    """

    dx = abs(x[1] - x[0])

    return 2.0 + (2 * (dx ** 2) * m / hbar ** 2) * V


def _M(x, V, m, hbar):
    """
    This function calculates the Matrix solution.

    Parameters:
    -----------
    x: array of size N
       linearly spaced points defining spatial domain of well (m)
    V: array of size N
       values of potential for each corresponding x value (J)
    m: array of size N
       electron mass (kg)
    hbar: float
          Planck's constant / 2 pi (Js)

    Out:
    ----
    M: array of size NxN
       array of calculated M values. Eigenvalues of M are energy
       eigenfunctions of potential. Eigenvectors of M are corresponding
       wavefunctions.
    """

    N = len(x)
    dx = abs(x[1] - x[0])

    beta = _beta(x, V, m, hbar)
    M = np.diagflat(beta) - np.eye(N, k=1) - np.eye(N, k=-1)

    return M * hbar ** 2 / (2.0 * m * dx ** 2)


def energy_eigenvalues(x, V, m=None, hbar=constants.hbar):
    """
    This function calculates the energy eigenvalues of a given potential
    over a domain of x.

    Parameters:
    -----------
    x: array of size N
       linearly spaced points defining spatial domain of well (m)
    V: array of size N
       values of potential for each corresponding x value (J)
    m: array of size N
       electron mass (kg)
    hbar: float
          Planck's constant / 2 pi (Js)

    Out:
    ----
    eigenvalues: array of size N
                 array of N energy eigenvalues for given potential.
                 Energy eigenvalues are in Joules.
    """

    if len(x) < 2:
        raise ValueError('Spatial coordinates must be more than a '
                         'single point')
    if len(x) != len(V):
        raise ValueError('array x and array V must have same length')
    if m is None:
        m = np.ones_like(x) * constants.m_e
    else:
        if len(x) != len(m):
            raise ValueError('array x and array m must have same '
                             'length')

    M = _M(x, V, m, hbar)
    eigenvalues = LA.eigvalsh(M, UPLO='U')

    return eigenvalues


def wavefunctions(x, V, m=None, hbar=constants.hbar):
    """
    This function calculates the wavefunctions for a given potential
    over a domain of x.

    Parameters:
    -----------
    x: array of size N
       linearly spaced points defining spatial domain of well (m)
    V: array of size N
       values of potential for each corresponding x value (J)
    m: float
       electron mass (kg)
    hbar: float
          Planck's constant / 2 pi (Js)

    Out:
    ----
    eigenvectors: array of size NxN
                  array of N wavefunctions. eigenvectors[:,0] is first
                  eigenvector along x, eigenvectors[:1] is 2nd etc.
    """

    if len(x) < 2:
        raise ValueError('Spatial coordinates must be more than a '
                         'single point')
    if len(x) != len(V):
        raise ValueError('array x and array V must have same length')
    if m is None:
        m = np.ones_like(x) * constants.m_e
    else:
        if len(x) != len(m):
            raise ValueError('array x and array m must have same '
                             'length')

    M = _M(x, V, m, hbar)
    _, eigenvectors = LA.eigh(M, UPLO='U')

    # Uncomment if eigvectors need flipping
    # for i in range(len(x)):
    #     if eigenvectors[0, i] > eigenvectors[1, i]:
    #         eigenvectors[:, i] *= -1

    return eigenvectors


def solve_schrodinger(x, V, m=None, hbar=constants.hbar):
    """
    This function calculates the energy eigenvalues and corresponding
    wavefunctions for a given potential over a domain of x.

    Parameters:
    -----------
    x: array of size N
       linearly spaced points defining spatial domain of well (m)
    V: array of size N
       values of potential for each corresponding x value (J)
    m: float
       electron mass (kg)
    hbar: float
          Planck's constant / 2 pi (Js)

    Out:
    ----
    eigenvalues: array of size N
                 array of N energy eigenvalues for given potential.
                 Energy eigenvalues are in Joules.
    eigenvectors: array of size NxN
                  array of N wavefunctions. eigenvectors[:,0] is first
                  eigenvector along x, eigenvectors[:1] is 2nd etc.
    """

    if len(x) < 2:
        raise ValueError('Spatial coordinates must be more than a '
                         'single point')
    if len(x) != len(V):
        raise ValueError('array x and array V must have same length')
    if m is None:
        m = np.ones_like(x) * constants.m_e
    else:
        if len(x) != len(m):
            raise ValueError('array x and array m must have same '
                             'length')

    M = _M(x, V, m, hbar)
    eigenvalues, eigenvectors = LA.eigh(M, UPLO='U')

    # Uncomment if eigvectors need flipping
    # for i in range(len(x)):
    #     if eigenvectors[0, i] > eigenvectors[1, i]:
    #         eigenvectors[:, i] *= -1

    return eigenvalues, eigenvectors


if __name__ == '__main__':
    import matplotlib.pyplot as plt

    """DEGENERATE STATES PRESENT"""
    # x = np.linspace(-1e-8, 1e-8, 1000)
    # V = np.ones_like(x) * 1.6e-19 * 200e-3
    # V[abs(x) < 6e-9] = 0
    # V[abs(x) < 2e-9] = 1.6e-19 * 100e-3
    #
    # eigval, eigvect = solve_schrodinger(x, V)
    # plt.figure()
    # plt.plot(x/1e-9, ((eigvect[:,0:5])*0.5e-19 + eigval[:5])/constants.eV)
    # plt.plot(x/1e-9, V/constants.eV, 'k')
    # plt.xlabel('z (nm)')
    # plt.ylabel('E (eV)')
    # plt.show()

    """QCL"""
    # x = np.linspace(-1e-8, 1e-8, 1000)
    # V = np.ones_like(x) * 1.6e-19 * 200e-3
    # V[100:150] = 0
    # V[200:300] = 0
    # V[375:400] = 0
    # V[450:500] = 0
    # V[525:575] = 0
    # V[650:700] = 0
    # V[825:900] = 0
    # V[925:975] = 0
    #
    # V -= 0.5*x*1e-11
    # V -= min(V)
    #
    # eigval, eigvect = solve_schrodinger(x, V)
    # plt.figure()
    # plt.plot(x, (eigvect[:, 0:10]) * 0.5e-19 + eigval[:10])
    # plt.plot(x, V, 'k')
    # plt.show()

    """OVERLAP INTEGRAL - CONST MASS VS VARYING MASS"""
    import materials
    from anderson import anderson
    x = np.linspace(0, 30e-9, 1000)
    mat = [materials.AlGaAs(0.2, 10), materials.GaAs(10),
           materials.AlGaAs(0.2, 10)]
    V = anderson(x, mat, np.array([1e-8, 1e-8, 1e-8]))

    m_const = np.ones_like(x) * (0.063) * constants.m_e #+ 0.083*0.2) * constants.m_e
    m_change = np.ones_like(x) * (0.063 + 0.083*0.2)
    m_change[x >= 10e-9] = 0.063
    m_change[x >= 20e-9] = (0.063 + 0.083*0.2)
    m_change *= constants.m_e

    eigval_const, eigvect_const = solve_schrodinger(x, V, m_const)
    eigval_change, eigvect_change = solve_schrodinger(x, V, m_change)

    plt.figure()
    plt.plot(x/1e-9, V/constants.eV, 'k')
    plt.plot(x / 1e-9, (eigvect_const[:, 0] + eigval_const[0] / constants.eV),
             label='const mass')
    plt.plot(x / 1e-9, (eigvect_const[:, 1] + eigval_const[1] / constants.eV))
    plt.plot(x / 1e-9, (eigvect_const[:, 2] + eigval_const[2] / constants.eV))
    plt.gca().set_prop_cycle(None)
    plt.plot(x / 1e-9, eigvect_change[:, 1] + eigval_change[1] / constants.eV,
             '--', label='changing mass')
    plt.plot(x / 1e-9, eigvect_change[:, 2] + eigval_change[2] / constants.eV, '--')
    plt.plot(x / 1e-9, eigvect_change[:, 3] + eigval_change[3] / constants.eV, '--')

    plt.legend(loc=1)

    plt.xlabel('z (nm)')
    plt.ylabel('E (eV)')
    # plt.show()

    for i in range(4):
        top = np.trapz(eigvect_const[:, i]*eigvect_change[:, i+1], x)
        bottom = np.trapz(eigvect_const[:, i]**2, x)
        print(f'Integral overlap {i+1} = {top/bottom:.3g}')

    print('\nConst\tChange\tdE  \tchange/const')
    for i in range(3):
        print(f'{eigval_const[i]/constants.eV:.3f}\t'
              f'{eigval_change[i+1]/constants.eV:.3f}\t'
              f'{abs(eigval_const[i] - eigval_change[i+1])/constants.eV:.3f}\t'
              f'{eigval_change[i+1]/eigval_const[i]:.3f}')


    # Determine x**2 dependence? Straight line if direct x**2 dependence
    # plt.figure()
    # plt.loglog(eigval_const, label='Constant mass')
    # plt.loglog(eigval_change[1:], label='Changing mass') # ignore E[0]
    # plt.ylabel('log(E)')
    # plt.title('All E states')
    # plt.grid()
    # plt.legend(loc='best')
    # States in well - only 3 points so kind of useless
    # plt.figure()
    # plt.loglog(eigval_const[:4], label='Constant mass')
    # plt.loglog(eigval_change[1:5], label='Changing mass')
    # plt.ylabel('log(E)')
    # plt.title('E states < max(V)')
    # plt.grid()
    # plt.legend(loc='best')

    """No wavefunc in well test - changing mass has weird effects"""
    import materials
    from anderson import anderson

    x = np.linspace(0, 40e-9, 1000)
    mat = [materials.AlGaAs(0.2, 10), materials.GaAs(10),
           materials.InAs(10), materials.AlGaAs(0.2, 10)]
    V = anderson(x, mat, np.array([1e-8, 1e-8, 1e-8, 1e-8]))

    m_const = np.ones_like(x) * 0.023 * constants.m_e
    # m_const = np.ones_like(x) *(0.063 + 0.083 *0.2) * constants.m_e
    m_change = np.ones_like(x) * (0.063 + 0.083 * 0.2)
    m_change[x >= 10e-9] = 0.063
    m_change[x >= 20e-9] = 0.023
    m_change[x >= 30e-9] = (0.063 + 0.083 * 0.2)
    m_change *= constants.m_e

    eigval_const, eigvect_const = solve_schrodinger(x, V, m_const)
    eigval_change, eigvect_change = solve_schrodinger(x, V, m_change)

    plt.figure()
    plt.plot(x/1e-9, V/constants.eV, 'k')
    plt.plot(x/1e-9, eigvect_const[:,0] + eigval_const[0]/constants.eV, label='const')
    plt.plot(x/1e-9, eigvect_const[:,1] + eigval_const[1]/constants.eV)
    plt.plot(x/1e-9, eigvect_const[:,2] + eigval_const[2]/constants.eV)
    plt.plot(x/1e-9, eigvect_const[:,3] + eigval_const[3]/constants.eV)
    # plt.plot(x/1e-9, eigvect_const[:,4] + eigval_const[4]/constants.eV)
    # plt.plot(x/1e-9, eigvect_const[:,5] + eigval_const[5]/constants.eV)
    # plt.plot(x/1e-9, eigvect_const[:,6] + eigval_const[6]/constants.eV)
    plt.gca().set_prop_cycle(None)

    plt.plot(x/1e-9, eigvect_change[:,2] + eigval_change[
        2]/constants.eV,
             '--', label='change')
    plt.plot(x / 1e-9, eigvect_change[:,3] + eigval_change[
        3]/constants.eV,
             '--')
    plt.plot(x / 1e-9, eigvect_change[:,4] + eigval_change[
        4]/constants.eV,
             '--')
    plt.plot(x / 1e-9, eigvect_change[:,5] + eigval_change[
        5]/constants.eV,
             '--')
    plt.legend(loc='best')


    plt.show()
