from numpy import pi


class ParametryMaterialowe:
    """
    Klasa, w której zawarte są parametry materiałowe.
    """

    mu0H0 = 0.1
    gamma = 176.6e9
    H0 = mu0H0 / (4e-7 * pi)
    mu0 = 4e-7 * pi
    d = 30e-9
    a = 16*90e-9
    b = 16*90e-9
    x = 0

    MoFe = 1.752e6
    MoNi = 0.484e6
    lNi = 7.64e-9 * 7.64e-9
    lFe = 3.3e-9 * 3.3e-9
    r = 0
    MoCo = 1.445e6
    MoPy = 0.658e6
    ACo = 3e-11
    APy = 1.1e-11

    MoA = MoCo # rdzeń
    MoB = MoPy # wypełnienie

    #lA = lNi
    #lB = lFe

    AA = APy
    AB = ACo

    """
    mu0 = 4e-7 * pi
    H0 = 0.1 / mu0
    mu0H0 = H0 * mu0
    gamma = 176e9
    d = 20e-9
    a = 400e-9
    b = 400e-9
    MoA = 1.752e6
    MoB = 0.484e6
    r = 160e-9
    x = 0
    lA = 3.3e-9 * 3.3e-9
    lB = 7.64e-9 * 7.64e-9
    """
    """
    mu0H0 = 0.074
    gamma = 176e9
    H0 = mu0H0 / (4e-7 * pi)
    mu0 = 4.e-7 * pi
    d = 10e-9
    a = 960e-9
    b = 960e-9
    MoPy = 0.8e6
    MoCo = 0
    APy = 1.1e-11
    ACo = 0
    x = 0.
    """
    def __init__(self, ilosc_wektorow, typ_pola_wymiany):
        self.ilosc_wektorow = ilosc_wektorow
        self.typ_pola_wymiany = typ_pola_wymiany

        if typ_pola_wymiany == 'I':
            self.lA = 2 * self.ACo / (self.mu0 * self.MoCo)
            self.lB = 2. * self.APy / (self.mu0 * self.MoPy)
        else:
            self.lA = 2. * self.ACo / (self.mu0 * self.MoCo ** 2)
            self.lB = 2. * self.APy / (self.mu0 * self.MoPy ** 2)

        """
        if typ_pola_wymiany == 'I':
            self.lCo = 0
            self.lPy = 2. * self.APy / (self.mu0 * self.MoPy)
        else:
            self.lCo = 0
            self.lPy = 2. * self.APy / (self.mu0 * self.MoPy ** 2)
        """