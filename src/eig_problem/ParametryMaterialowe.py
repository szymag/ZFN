from numpy import pi


class ParametryMaterialowe:
    """
    Klasa, w której zawarte są parametry materiałowe.
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

    def __init__(self, ilosc_wektorow, typ_pola_wymiany):
        self.ilosc_wektorow = ilosc_wektorow
        self.typ_pola_wymiany = typ_pola_wymiany
        if typ_pola_wymiany == 'I':
            self.lCo = 2 * self.ACo / (self.mu0 * self.MoCo)
            self.lPy = 2 * self.APy / (self.mu0 * self.MoPy)
        else:
            self.lCo = 2 * self.ACo / (self.mu0 * self.MoCo)
            self.lPy = 2 * self.APy / (self.mu0 * self.MoPy)
