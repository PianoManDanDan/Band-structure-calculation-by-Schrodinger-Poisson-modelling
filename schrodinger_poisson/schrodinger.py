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
    eigenvalues = LA.eigvals(M)

    return eigenvalues.sort()


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
    eigenvalues, eigenvectors = LA.eig(M)

    # sort eigenvectors into ascending order
    idx = eigenvalues.argsort()
    eigenvectors = eigenvectors[:, idx]

    # Uncomment if eigvectors need flipping
    # for i in range(len(x)):
    #     if eigenvectors[0, i] > eigenvectors[1, i]:
    #         eigenvectors[:, i] *= -1

    return eigenvectors




def myMatrix1(x, V, m, hbar=constants.hbar): ## Doesn't work

    dx = abs(x[1] - x[0])

    Ma = 1 / (dx**2 * m[:-1])
    Mc = 1 / (dx**2 * np.append(m, m[-1])[2:])
    Md_top = - hbar**2 * m - hbar**2 * np.append(m, m[-1])[1:] \
             - 2 * dx**2 * m * np.append(m, m[-1])[1:] * V
    Md_bottom = hbar**2 * dx**2 * m * np.append(m, m[-1])[1:]
    Md = Md_top / Md_bottom

    M = np.diagflat(Md) + np.diagflat(Mc, 1) + np.diagflat(Ma, -1)

    return M * -hbar**2 / 2

def myMatrix2(x, V, m, hbar=constants.hbar):

    dx = abs(x[1] - x[0])

    Ma = 1 / (dx**2 * m[1:])
    Mc = 1 / (dx**2 * m[1:])
    Md_top = - hbar**2 * m - hbar**2 * np.append(m, m[-1])[1:] - 2 * \
             dx**2 * m * np.append(m, m[-1])[1:] * V
    Md_bottom = hbar**2 * dx**2 * m * np.append(m, m[-1])[1:]
    Md = Md_top / Md_bottom

    M = np.diagflat(Md) + np.diagflat(Mc, 1) + np.diagflat(Ma, -1)

    return M * -hbar**2 / 2

def myMatrix3(x, V, m, hbar=constants.hbar): ## Doesn't work

    dx = abs(x[1] - x[0])

    Ma = 1 / (dx**2 * m[:-1])
    Mc = 1 / (dx**2 * m[1:])
    Md_top = - hbar**2 * m - hbar**2 * np.append(m, m[-1])[1:] - 2 * \
             dx**2 * m * np.append(m, m[-1])[1:] * V
    Md_bottom = hbar**2 * dx**2 * m * np.append(m, m[-1])[1:]
    Md = Md_top / Md_bottom

    M = np.diagflat(Md) + np.diagflat(Mc, 1) + np.diagflat(Ma, -1)

    return M * -hbar**2 / 2

def myMatrix4(x, V, m, hbar=constants.hbar): ## Doesn't work

    dx = abs(x[1] - x[0])

    Ma = 1 / (dx**2 * m[1:])
    Mc = 1 / (dx**2 * np.append(m, m[-1])[2:])
    Md_top = - hbar**2 * m - hbar**2 * np.append(m, m[-1])[1:] \
             - 2 * dx**2 * m * np.append(m, m[-1])[1:] * V
    Md_bottom = hbar**2 * dx**2 * m * np.append(m, m[-1])[1:]
    Md = Md_top / Md_bottom

    M = np.diagflat(Md) + np.diagflat(Mc, 1) + np.diagflat(Ma, -1)

    return M * -hbar**2 / 2

def TristanMatrix(x, V, m, hbar=constants.hbar):
    dx = abs(x[1] - x[0])

    Ma = 2*np.append(m, [m[-1]])[2:]*m[:-1]
    Ma -= m[1:]*(m[:-1] - np.append(m, [m[-1]])[2:])*dx

    Mc = 2*m[1:]*np.append([m[0]], m)[:-2]
    Mc += m[:-1]*(np.append([m[0]], m)[:-2] - m[1:])*dx

    Md = (-4*dx**2 / hbar**2)
    Md *= (np.append(m, [m[-1]])[1:] * m * np.append([m[0]], m)[:-1])*V
    Md -= 4*(np.append(m, [m[-1]])[1:] * np.append([m[0]], m)[:-1])

    abottom = np.append(m, [m[-1]])[2:]*m[:-1] * m[1:]
    dbottom = np.append(m, [m[-1]])[1:] * m * np.append([m[0]], m)[:-1]
    cbottom = m[1:]*np.append([m[0]], m)[:-2] * m[:-1]

    Ma /= abottom
    Md /= dbottom
    Mc /= cbottom

    M = np.diagflat(Md) + np.diagflat(Mc, 1) + np.diagflat(Ma, -1)

    return M * -hbar**2 / (4*dx**2)




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
    # M = myMatrix1(x, V, m, hbar)
    # M = myMatrix2(x, V, m, hbar)
    # M = myMatrix3(x, V, m, hbar)
    # M = myMatrix4(x, V, m, hbar)
    # M = TristanMatrix(x, V, m, hbar)
    eigenvalues, eigenvectors = LA.eigh(M)

    # sort eigenvalues and eigenvectors into ascending order
    idx = eigenvalues.argsort()
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]

    # comment if eigvectors don't need flipping
    for i in range(len(x)):
        if eigenvectors[0, i] > eigenvectors[1, i]:
            eigenvectors[:, i] *= -1

    return eigenvalues.real, eigenvectors.real






if __name__ == '__main__':
    import matplotlib.pyplot as plt

    """DEGENERATE STATES PRESENT"""
    x = np.linspace(-1e-8, 1e-8, 1000)
    V = np.ones_like(x) * 1.6e-19 * 200e-3
    V[abs(x) < 6e-9] = 0
    V[abs(x) < 2e-9] = 1.6e-19 * 100e-3

    eigval, eigvect = solve_schrodinger(x, V)
    print(eigval[:5] *1000/ constants.eV)
    plt.figure()
    plt.plot(x/1e-9, ((eigvect[:,0:5])*0.5e-19 + eigval[
                                                 :5])*1000/constants.eV)
    plt.plot(x/1e-9, V*1000/constants.eV, 'k')
    plt.xlabel('z (nm)')
    plt.ylabel('E (meV)')
    plt.show()

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
    # import materials
    # from anderson import anderson
    # x = np.linspace(0, 50e-9, 1000)
    # mat = [materials.AlGaAs(0.5, 10), materials.GaAs(10),
    #        materials.AlGaAs(0.5, 10)]
    # V = anderson(x, mat, np.array([2e-8, 1e-8, 2e-8]))
    #
    # m_const = np.ones_like(x) * (0.063) * constants.m_e #+ 0.083*0.2) * constants.m_e
    # m_change = np.ones_like(x) * (0.063 + 0.083*0.2)
    # m_change[x >= 10e-9] = 0.063
    # m_change[x >= 20e-9] = (0.063 + 0.083*0.2)
    # m_change *= constants.m_e
    #
    # eigval_const, eigvect_const = solve_schrodinger(x, V, m_const)
    # eigval_change, eigvect_change = solve_schrodinger(x, V, m_change)
    #
    # plt.figure()
    # plt.plot(x/1e-9, V/constants.eV, 'k')
    # plt.plot(x / 1e-9, (eigvect_const[:, 0] + eigval_const[0] / constants.eV),
    #          label='const mass')
    # plt.plot(x / 1e-9, (eigvect_const[:, 1] + eigval_const[1] / constants.eV))
    # plt.plot(x / 1e-9, (eigvect_const[:, 2] + eigval_const[2] / constants.eV))
    # plt.gca().set_prop_cycle(None)
    # plt.plot(x / 1e-9, eigvect_change[:, 0] + eigval_change[0] /
    #          constants.eV,
    #          '--', label='changing mass')
    # plt.plot(x / 1e-9, eigvect_change[:, 1] + eigval_change[1] /
    #          constants.eV, '--')
    # plt.plot(x / 1e-9, eigvect_change[:, 2] + eigval_change[2] /
    #          constants.eV, '--')
    #
    # plt.legend(loc=1)
    #
    # plt.xlabel('z (nm)')
    # plt.ylabel('E (eV)')
    # # plt.show()
    #
    # for i in range(4):
    #     top = np.trapz(eigvect_const[:, i]*eigvect_change[:, i], x)
    #     bottom = np.trapz(eigvect_const[:, i]**2, x)
    #     print(f'Integral overlap {i+1} = {top/bottom:.3g}')
    #
    # print('\nConst\tChange\tdE  \tchange/const')
    # for i in range(3):
    #     print(f'{eigval_const[i]/constants.eV:.3f}\t'
    #           f'{eigval_change[i+1]/constants.eV:.3f}\t'
    #           f'{abs(eigval_const[i] - eigval_change[i+1])/constants.eV:.3f}\t'
    #           f'{eigval_change[i+1]/eigval_const[i]:.3f}')

    """Big discontinuities test"""
    # import materials
    # from anderson import anderson
    #
    # x = np.linspace(0, 40e-9, 1000)
    # mat = [materials.AlGaAs(0.2, 10), materials.GaAs(10),
    #        materials.InAs(10), materials.AlGaAs(0.2, 10)]
    # V = anderson(x, mat, np.array([1e-8, 1e-8, 1e-8, 1e-8]))
    #
    # m_const = np.ones_like(x) * 0.023 * constants.m_e
    # # m_const = np.ones_like(x) *(0.063 + 0.083 *0.2) * constants.m_e
    # m_change = np.ones_like(x) * (0.063 + 0.083 * 0.2)
    # m_change[x >= 10e-9] = 0.063
    # m_change[x >= 20e-9] = 0.023
    # m_change[x >= 30e-9] = (0.063 + 0.083 * 0.2)
    # m_change *= constants.m_e
    #
    # eigval_const, eigvect_const = solve_schrodinger(x, V, m_const)
    # eigval_change, eigvect_change = solve_schrodinger(x, V, m_change)
    #
    # plt.figure()
    # plt.plot(x/1e-9, V/constants.eV, 'k')
    # plt.plot(x/1e-9, eigvect_const[:,0] + eigval_const[0]/constants.eV, label='const')
    # plt.plot(x/1e-9, eigvect_const[:,1] + eigval_const[1]/constants.eV)
    # plt.plot(x/1e-9, eigvect_const[:,2] + eigval_const[2]/constants.eV)
    # plt.plot(x/1e-9, eigvect_const[:,3] + eigval_const[3]/constants.eV)
    # # plt.plot(x/1e-9, eigvect_const[:,4] + eigval_const[4]/constants.eV)
    # # plt.plot(x/1e-9, eigvect_const[:,5] + eigval_const[5]/constants.eV)
    # # plt.plot(x/1e-9, eigvect_const[:,6] + eigval_const[6]/constants.eV)
    # plt.gca().set_prop_cycle(None)
    #
    # plt.plot(x/1e-9, eigvect_change[:,0] + eigval_change[
    #     0]/constants.eV,
    #          '--', label='change')
    # plt.plot(x / 1e-9, eigvect_change[:,1] + eigval_change[
    #     1]/constants.eV,
    #          '--')
    # plt.plot(x / 1e-9, eigvect_change[:,2] + eigval_change[
    #     2]/constants.eV,
    #          '--')
    # plt.plot(x / 1e-9, eigvect_change[:,3] + eigval_change[
    #     3]/constants.eV,
    #          '--')
    # plt.legend(loc='best')

    """Small changes in mass"""
    # x = np.linspace(0, 30e-9, 100)
    # mat = [materials.AlGaAs(0.2, 10), materials.GaAs(10),
    #        materials.AlGaAs(0.2, 10)]
    # V = anderson(x, mat, np.array([1e-8, 1e-8, 1e-8]))
    #
    # m_const = np.ones_like(x) * (0.063) * constants.m_e
    # m_change = np.ones_like(x) * (0.063 + 0.01)
    # m_change[x >= 10e-9] = 0.063
    # m_change[x >= 20e-9] = (0.063 + 0.01)
    # m_change *= constants.m_e
    #
    # eigval_const, eigvect_const = solve_schrodinger(x, V, m_const)
    # eigval_change, eigvect_change = solve_schrodinger(x, V, m_change)
    #
    # plt.figure()
    # plt.plot(x / 1e-9, V / constants.eV, 'k')
    # plt.plot(x / 1e-9,
    #          (eigvect_const[:, 0] + eigval_const[0] / constants.eV),
    #          label='const mass')
    # plt.plot(x / 1e-9,
    #          (eigvect_const[:, 1] + eigval_const[1] / constants.eV))
    # plt.plot(x / 1e-9,
    #          (eigvect_const[:, 2] + eigval_const[2] / constants.eV))
    # plt.gca().set_prop_cycle(None)
    # plt.plot(x / 1e-9,
    #          eigvect_change[:, 0] + eigval_change[0] / constants.eV,
    #          '--', label='changing mass')
    # plt.plot(x / 1e-9,
    #          eigvect_change[:, 1] + eigval_change[1] / constants.eV,
    #          '--')
    # plt.plot(x / 1e-9,
    #          eigvect_change[:, 2] + eigval_change[2] / constants.eV,
    #          '--')
    # plt.title('small dV')
    #
    # plt.legend(loc=1)
    # plt.xlabel('z (nm)')
    # plt.ylabel('E (eV)')

    """Cunha, Christansen 2017 test"""
    # z = np.linspace(-np.pi/2, np.pi/2, 1000)
    # # z = np.linspace(-1.2, 1.2, 1000)
    # x = np.arccosh(1/np.cos(z))
    # # x = np.linspace(-0.5, 0.5, 1000)
    # Vx = np.zeros_like(x)
    # Vz = 0.5 + 0.75*np.tan(z)**2
    # mx = constants.m_e / np.cosh(x)**2
    # mz = constants.m_e * np.cos(z)
    #
    #
    # for n in range(1, 5):
    #     analyticE = (constants.hbar / (2 * constants.m_e)) * (n * (n+1))
    #     print(f'E_{n} = {analyticE * constants.eV:.3e} J')
    # print('\n')
    #
    #
    # evals, evects = solve_schrodinger(x, Vx, mx)
    #
    # plt.figure()
    # plt.plot(x[x<0.5], Vx[x<0.5], 'k')
    # for i in range(3):
    #     plt.plot(x[x<0.5], evals[i] + evects[x<0.5,i])
    #     print(f'Calc E_{i+1} = {evals[i]:.3e} eV')
    # plt.grid()



    # z = np.linspace(-np.pi/2, np.pi/2, 1000)
    # # z = np.linspace(-1.2, 1.2, 1000)
    # x = np.arccosh(1/np.cos(z))
    # # x = np.linspace(-0.5, 0.5, 1000)
    # Vx = - (3 * constants.hbar**2 / (8 * constants.m_e)) * np.sinh(
    #     x)**2 - constants.hbar**2 / (4 * constants.m_e)
    # Vz = - (3 * constants.hbar**2 / (8 * constants.m_e)) * np.tan(
    #     z)**2 - constants.hbar**2 / (4 * constants.m_e)
    # mx = constants.m_e / np.cosh(x)**2
    # mz = constants.m_e * np.cos(z)
    #
    # for n in range(1, 5):
    #     analyticE = (constants.hbar**2 * np.pi**2 / (8 * constants.m_e)) * n**2
    #     print(f'E_{n} = {analyticE * constants.eV:.3e} J')
    # print('\n')
    #
    #
    # evals, evects = solve_schrodinger(z, Vz, mz)
    #
    # plt.figure()
    # plt.plot(z, Vz, 'k')
    # for i in range(3):
    #     plt.plot(z, evals[i] + evects[:,i])
    #     print(f'Calc E_{i+1} = {evals[i]:.3e} eV')
    # plt.grid()

    """Analytic vs matrix"""
    # def E_n(n, L):
    #     top = n**2 * np.pi**2 * constants.hbar**2
    #     bottom = 2 * constants.m_e *L**2
    #     return top / bottom
    #
    # L = 10e-9
    #
    # analytic = np.zeros(5)
    # for i in range(5):
    #     analytic[i] = E_n(i+1, L)
    # analytic /= (constants.eV * 1e-3)
    #
    # x = np.linspace(0, L, 1000)
    # V = np.zeros_like(x)
    # m = np.ones_like(x) * constants.m_e
    # evals, evects = solve_schrodinger(x, V, m)
    # evals /= (constants.eV * 1e-3)
    #
    # print(f'Anal\tMatr\tdE')
    # for i in range(5):
    #     print(f'{analytic[i]:.3f}\t{evals[i]:.3f}\t'
    #           f'{abs(analytic[i] - evals[i]):.3f}')
    #
    # plt.figure()
    # plt.plot(x, V, 'k')
    # for i in range(5):
    #     plt.plot(x, evects[:,i]*1.2e2 + evals[i])

    """Parabolic Well"""
    def E_n(n, omega):
        return (n+0.5) * constants.hbar * omega

    omega = np.sqrt(1 / constants.m_e)

    analytic = np.zeros(5)
    for i in range(5):
        analytic[i] = E_n(i, omega)
    analytic /= constants.eV

    x = np.linspace(0, 3e-9, 1000)
    V = 0.5 * omega**2 * constants.m_e * (x-max(x)/2)**2
    m = np.ones_like(x) * constants.m_e
    evals, evects = solve_schrodinger(x, V, m)
    evals /= constants.eV

    print(f'Anal\tMatr\tdE')
    for i in range(5):
        print(f'{analytic[i]:.3f}\t{evals[i]:.3f}\t'
              f'{(analytic[i] - evals[i])*1000:.3f}')

    plt.figure()

    plt.plot(x[V<3.5*constants.eV]/1e-9,
             V[V<3.5*constants.eV]/constants.eV, 'k')
    for i in range(5):
        plt.plot(x/1e-9, evects[:,i] + evals[i])

    plt.grid()
    plt.xlabel('x (nm)')
    plt.ylabel('E (eV)')
    plt.show()
