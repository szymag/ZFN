from numpy import pi


class ParametryMaterialowe:
    """
    Klasa, w której zawarte są parametry materiałowe.
    """
    """
    mu0H0 = 0.2
    gamma = 194.6e9
    H0 = mu0H0 / (4e-7 * pi)
    mu0 = 4e-7 * pi
    d = 4e-9
    a = 30e-9
    b = 30e-9
    MoCo = 1.15e6
    MoPy = 0.658e6
    ACo = 2.88e-11
    APy = 1.1e-11
    r = 14e-9
    x = 0

    """
    mu0 = 4e-7 * pi
    H0 = 0.1 / mu0
    mu0H0 = H0 * mu0
    gamma = 176e9
    d = 20e-9
    a = 400e-9
    b = 400e-9
    MoCo = 1.752e6
    MoPy = 0.484e6
    r = 160e-9
    x = 0
    lCo = 3.3e-9 * 3.3e-9
    lPy = 7.64e-9 * 7.64e-9

    """
    mu0H0 = 0.5
    gamma = 194.6e9
    H0 = mu0H0 / (4e-7 * pi)
    mu0 = 4.e-7 * pi
    d = 10e-9
    a = 30e-9
    b = 30e-9
    MoPy = 0.658e6
    MoCo = 0
    APy = 1.1e-11
    ACo = 0
    x = 0.
"""
    def __init__(self, ilosc_wektorow, typ_pola_wymiany):
        self.ilosc_wektorow = ilosc_wektorow
        self.typ_pola_wymiany = typ_pola_wymiany
"""
        if typ_pola_wymiany == 'I':
            self.lCo = 2 * self.ACo / (self.mu0 * self.MoCo)
            self.lPy = 2. * self.APy / (self.mu0 * self.MoPy)
        else:
            self.lCo = 2. * self.ACo / (self.mu0 * self.MoCo ** 2)
            self.lPy = 2. * self.APy / (self.mu0 * self.MoPy ** 2)
        """
"""
        if typ_pola_wymiany == 'I':
            self.lCo = 0
            self.lPy = 2. * self.APy / (self.mu0 * self.MoPy)
        else:
            self.lCo = 0
            self.lPy = 2. * self.APy / (self.mu0 * self.MoPy ** 2)
        """