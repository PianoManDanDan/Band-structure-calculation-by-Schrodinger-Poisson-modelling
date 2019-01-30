"""
This is a module for solving Schrodinger's equation using the matrix
method. 

The main function of the module (solve_schrodinger) returns the energy
eigenvalues and corresponding wavefunctions for a potential, V, in a
domain, x.
"""

__version__ = '0.1'
__author__ = 'Daniel Martin'
__all__ = ['_beta', '_M', 'energy_eigenvalues',
        'wavefunctions', 'solve_schrodinger']

import numpy as np
from numpy import linalg as LA
from scipy import constants


def _beta(x, V, m, hbar):
    """
    This function calculates the beta value matrix.

    Parameters:
    -----------
    x: array of size N
       linearly spaced, spatial domain of well (m)
    V: array of size N
       values of potential for each corresponding x value (J)
    m: float
       electron mass (kg)
    hbar: float
          Planck's constant / 2 pi (Js)

    Out:
    ----
    beta: array of size N
          array of beta values to go into calculating M
    """

    dx = abs(x[1]-x[0])

    return 2.0 + (2*(dx**2)*m / hbar**2)*V


def _M(x, V, m, hbar):
    """
    This function calculates the Matrix solution.

    Parameters:
    -----------
    x: array of size N
       linearly spaced, spatial domain of well (m)
    V: array of size N
       values of potential for each corresponding x value (J)
    m: float
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
    dx = abs(x[1]-x[0])

    beta = _beta(x, V, m, hbar)
    M = np.diagflat(beta) - np.eye(N, k=1) - np.eye(N, k=-1)

    return M * hbar**2 / (2.0 * m * dx**2)


def energy_eigenvalues(x, V, m=constants.m_e, hbar=constants.hbar):
    """
    This function calculates the energy eigenvalues of a given potential
    over a domain of x.

    Parameters:
    -----------
    x: array of size N
       linearly spaced, spatial domain of well (m)
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
    """

    # Check x and V have same length
    if len(x) != len(V):
        raise ValueError('array x and array V must have same length')

    M = _M(x, V, m, hbar)
    eigenvalues = LA.eigvalsh(M, UPLO='U')

    return eigenvalues


def wavefunctions(x, V, m=constants.m_e, hbar=constants.hbar):
    """
    This function calculates the wavefunctions for a given potential
    over a domain of x.

    Parameters:
    -----------
    x: array of size N
       linearly spaced, spatial domain of well (m)
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

    # Check x and V have same length
    if len(x) != len(V):
        raise ValueError('array x and array V must have same length')

    M = _M(x, V, m, hbar)
    _, eigenvectors = LA.eigh(M, UPLO='U')

    for i in range(len(x)):
        if eigenvectors[0,i] > eigenvectors[1,i]:
            eigenvectors[:,i] *= -1

    return eigenvectors


def solve_schrodinger(x, V, m=constants.m_e, hbar=constants.hbar):
    """
    This function calculates the energy eigenvalues and corresponding
    wavefunctions for a given potential over a domain of x.

    Parameters:
    -----------
    x: array of size N
       linearly spaced, spatial domain of well (m)
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

    # Check x and V have same length
    if len(x) != len(V):
        raise ValueError('array x and array V must have same length')

    M = _M(x, V, m, hbar)
    eigenvalues, eigenvectors = LA.eigh(M, UPLO='U')

    for i in range(len(x)):
        if eigenvectors[0,i] > eigenvectors[1,i]:
            eigenvectors[:,i] *= -1

    return eigenvalues, eigenvectors




if __name__ == '__main__':
    x = np.linspace(-1e-8, 1e-8, 1000)
    V = np.zeros_like(x) 

    eigval, eigvect = solve_schrodinger(x, V)
