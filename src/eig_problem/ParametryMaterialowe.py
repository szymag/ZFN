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
    MoCo = 1.15e6
    MoPy = 0.658e6
    ACo = 2.88e-11
    APy = 1.1e-11
    lCo = 2 * ACo / (mu0 * MoCo)
    lPy = 2 * APy / (mu0 * MoPy)
    r = 7e-9
    x = 0

    def __init__(self, rozmiar_macierzy_blok):
        self.rozmiar_macierzy_blok = rozmiar_macierzy_blok
