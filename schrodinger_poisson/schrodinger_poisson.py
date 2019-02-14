"""
This is a module for solving the Schrodinger-Poisson band structure
model using the python scripts found in the rest of the
schrodinger_poisson module. This contains the self consistency test,
as well as other functions needed to perform the full calculation.
"""

__version__ = '0.1'
__author__ = 'Daniel Martin'
__all__ = ['self_consistent', 'solve_schrodinger_poisson']

import numpy as np
from scipy import constants

import schrodinger
import poisson
import diffusion


def self_consistent():
    """
    Determines whether a self consistent solution to the
    Schrodinger-Poisson solver has been reached.

    Parameters:
    -----------


    Out:
    ----
    consistent: boolean
                Boolean value determining whether a self consistent
                solution has been reached by the Schrodinger-Poisson
                solver.
    """

    consistent = False

    return consistent


def solve_schrodinger_poisson(x, V, volume_density,
                              relative_permittivity, m=constants.m_e,
                              hbar=constants.hbar, max_tests=1000):
    """
    Solves Schrodinger's and Poisson's equation in turn until a self
    consistent solution is reached.

    Parameters:
    -----------


    Out:
    ----

    """

    consistent = False
    tests = 0

    while not consistent and tests < max_tests:
        eigval, eigvect = schrodinger.solve_schrodinger(x, V, m, hbar)
        V = poisson.solve_poisson(x, V, eigvect[:,0], volume_density,
                                  relative_permittivity, q)
        consistent = self_consistent()
        tests += 1


if __name__ == '__main__':
    pass
