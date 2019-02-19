"""
This python module sets out classes for calculating the different
properties of common semiconductor materials.
"""

__version__ = '0.1'
__author__ = 'Daniel Martin'
__all__ = ['AlGaAs', 'Si', 'materials']

from scipy import constants


class AlGaAs:
    """
    class object for AlGaAs, GaAs and AlAs.

    Properties:
    -----------
    dielectric_constant: float
                         dielectric constant of material (units?)
    me: float
        Effective electron mass (kg)
    mlh: float
         Light hole mass (kg)
    mhh: float
         Heavy hole mass (kg)
    Eg: float
        Band gap energy (eV)
    """
    def __init__(self, x, T):
        """
        Sets up properties for AlGaAs, AlAs and GaAs.

        Parameters:
        -----------
        x: float
           A float between 0 and 1 inclusive. Determines proportion
           of Aluminium present in semiconductor material.
        T: float
           A float greater than zero. Temperature of material in Kelvin.
        """

        try:
            float(x)
            float(T)
        except TypeError:
            raise TypeError('Amount of Al in AlGaAs and Temperature '
                            'must be a float')
        if (x < 0) or (x > 1):
            raise ValueError('Amount of Al in AlGaAs must be between '
                             '0 and 1')
        if T < 0:
            raise ValueError('Temperature must be at least 0 Kelvin')

        self.dielectric_const = (12.9 - 2.84*x) * constants.epsilon_0
        self.mlh = (0.082 + 0.068*x) * constants.m_e
        self.mhh = (0.51 + 0.25*x) * constants.m_e
        if x == 0:
            self.me = 0.063 * constants.m_e
            self.Eg = 1.519 - 5.405e-4*T**2 / (T+204)
        # elif x == 1:          ## CANNOT FIND Eg or me FOR AlAs
            # self.me =
            # self.Eg =
        elif x < 0.45:
            self.me = (0.063 + 0.083*x) * constants.m_e
            self.Eg = 1.424 + 1.247*x
        else:
            # self.m_e = ??? - Not on IOFFE
            self.Eg = 1.9 + 0.125*x + 0.143*x**2


class InSb:
    """
    class object for AlGaAs, GaAs and AlAs.

    Properties:
    -----------
    dielectric_constant: float
                         dielectric constant of material (units?)
    me: float
        Effective electron mass (kg)
    mlh: float
         Light hole mass (kg)
    mhh: float
         Heavy hole mass (kg)
    Eg: float
        Band gap energy (eV)
        """

    def __init__(self, T):
        """
        Sets up properties for Si.

        Parameters
        ==========
        T: float
           A float greater than zero. Temperature of material in Kelvin.
        """

        try:
            float(T)
        except TypeError:
            raise TypeError('Temperature must be a float')
        if T < 0:
            raise ValueError('Temperature must be at least 0 Kelvin')

        self.dielectric_constant = 16.8 * constants.epsilon_0
        self.me = 0.014 * constants.m_e
        self.mlh = 0.015 * constants.m_e
        self.mhh = 0.43 * constants.m_e
        self.Eg = 0.24 - 6e-4 * T**2 / (T + 500)


class Si:
    """
    class object for AlGaAs, GaAs and AlAs.

    Properties:
    -----------
    dielectric_constant: float
                         dielectric constant of material (units?)
    me_long: float
             Longitudinal effective electron mass (kg)
    me_trans: float
              Transverse effective electron mass (kg)
    mlh: float
         Light hole mass (kg)
    mhh: float
         Heavy hole mass (kg)
    Eg: float
        Band gap energy (eV)
    """
    def __init__(self, T):
        """
        Sets up properties for Si.

        Parameters
        ==========
        T: float
           A float greater than zero. Temperature of material in Kelvin.
        """

        try:
            float(T)
        except TypeError:
            raise TypeError('Temperature must be a float')
        if T < 0:
            raise ValueError('Temperature must be at least 0 Kelvin')

        self.dielectric_constant = 11.7 * constants.epsilon_0
        self.me_long = 0.98 * constants.m_e
        self.me_trans = 0.19 * constants.m_e
        self.mlh = 0.16 * constants.m_e
        self.mhh = 0.49 * constants.m_e
        self.Eg = 1.17 - 4.73e-4 * T**2 / (T + 636)


materials = {'AlGaAs': AlGaAs, 'AlAs': AlGaAs, 'GaAs': AlGaAs,
             'InSb': InSb, 'Si': Si}


if __name__ == '__main__':
    AlGaAs1 = materials['AlGaAs'](0.5, 20)
    AlAs1 = materials['AlAs'](1, 20)
    GaAs1 = materials['GaAs'](0, 20)

    print(AlGaAs1.mhh)
    print(AlAs1.mhh)
    print(GaAs1.mhh)
