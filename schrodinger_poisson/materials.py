"""
This python module sets out classes for calculating the different
properties of common semiconductor materials.
"""

__version__ = '0.1'
__author__ = 'Daniel Martin'
__all__ = ['AlGaAs', 'InSb', 'InP', 'InAs', 'GaSb', 'materials']


class AlGaAs:
    """
    class object for AlGaAs, GaAs and AlAs.

    Properties:
    -----------
    dielectric_constant: float
                         dielectric constant of material (units?)
    me: float
        Effective electron mass (m0)
    mlh: float
         Light hole mass (m0)
    mhh: float
         Heavy hole mass (m0)
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

        self.dielectric_constant = 12.9 - 2.84*x
        self.mlh = 0.082 + 0.068*x
        self.mhh = 0.51 + 0.25*x
        if x == 0:
            self.me = 0.063
            self.Eg = 1.519 - 5.405e-4*T**2 / (T+204)
        elif x < 0.45:
            self.me = 0.063 + 0.083*x
            self.Eg = 1.424 + 1.247*x
        else:
            self.m_e = 0.063 + 0.083*0.45
            self.Eg = 1.9 + 0.125*x + 0.143*x**2


class InSb:
    """
    class object for InSb.

    Properties:
    -----------
    dielectric_constant: float
                         dielectric constant of material (units?)
    me: float
        Effective electron mass (m0)
    mlh: float
         Light hole mass (m0)
    mhh: float
         Heavy hole mass (m0)
    Eg: float
        Band gap energy (eV)
    """

    def __init__(self, _, T):
        """
        Sets up properties for InSb.

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

        self.dielectric_constant = 16.8
        self.me = 0.014
        self.mlh = 0.015
        self.mhh = 0.43
        self.Eg = 0.24 - 6e-4 * T**2 / (T + 500)


class InP:
    """
    class object for InP.

    Properties:
    -----------
    dielectric_constant: float
                         dielectric constant of material (units?)
    me: float
        Effective electron mass (m0)
    mlh: float
         Light hole mass (m0)
    mhh: float
         Heavy hole mass (m0)
    Eg: float
        Band gap energy (eV)
    """

    def __init__(self, _, T):
        """
        Sets up properties for InP.

        Parameters:
        -----------
        T: float
           A float greater than zero. Temperature of material in Kelvin.
        """

        try:
            float(T)
        except TypeError:
            raise TypeError('Temperature must be a float')
        if T < 0:
            raise ValueError('Temperature must be at least 0 Kelvin')

        self.dielectric_constant = 12.5
        self.me = 0.08
        self.mlh = 0.089
        self.mhh = 0.6
        self.Eg = 1.421 - 4.9e-4 * T**2 / (T + 327)


class InAs:
    """
    class object for InAs.

    Properties:
    -----------
    dielectric_constant: float
                         dielectric constant of material (units?)
    me: float
        Effective electron mass (m0)
    mlh: float
         Light hole mass (m0)
    mhh: float
         Heavy hole mass (m0)
    Eg: float
        Band gap energy (eV)
    """

    def __init__(self, _, T):
        """
        Sets up properties for InAs.

        Parameters:
        -----------
        T: float
           A float greater than zero. Temperature of material in Kelvin.
        """

        try:
            float(T)
        except TypeError:
            raise TypeError('Temperature must be a float')
        if T < 0:
            raise ValueError('Temperature must be at least 0 Kelvin')

        self.dielectric_constant = 15.15
        self.me = 0.023
        self.mlh = 0.026
        self.mhh = 0.41
        self.Eg = 0.415 - 2.76e-4 * T**2 / (T + 83)


class GaSb:
    """
    class object for GaSb.

    Properties:
    -----------
    dielectric_constant: float
                         dielectric constant of material (units?)
    me: float
        Effective electron mass (m0)
    mlh: float
         Light hole mass (m0)
    mhh: float
         Heavy hole mass (m0)
    Eg: float
        Band gap energy (eV)
    """

    def __init__(self, _, T):
        """
        Sets up properties for GaSb.

        Parameters:
        -----------
        T: float
           A float greater than zero. Temperature of material in Kelvin.
        """

        try:
            float(T)
        except TypeError:
            raise TypeError('Temperature must be a float')
        if T < 0:
            raise ValueError('Temperature must be at least 0 Kelvin')

        self.dielectric_constant = 15.7
        self.me = 0.041
        self.mlh = 0.05
        self.mhh = 0.4
        self.Eg = 0.813 - 3.78e-4 * T**2 / (T + 94)


materials = {'AlGaAs': AlGaAs, 'AlAs': AlGaAs, 'GaAs': AlGaAs,
             'InSb': InSb, 'InP': InP, 'InAs': InAs, 'GaSb': GaSb}


if __name__ == '__main__':
    AlGaAs = materials['AlGaAs'](0.5, 20)
    AlAs = materials['AlAs'](1, 20)
    GaAs = materials['GaAs'](0, 20)

    print(AlGaAs.mhh)
    print(AlAs.mhh)
    print(GaAs.mhh)
