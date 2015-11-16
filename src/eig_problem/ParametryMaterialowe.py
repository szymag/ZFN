from numpy import pi


class ParametryMaterialowe:

    mu0H0 = 0.001
    gamma = 194.6e9
    H0 = mu0H0 / (4e-7 * pi)
    d = 4e-9
    a = 30e-9
    MoCo = 1.15e6
    MoPy = 0.658e6
    ACo = 2.88e-11
    APy = 1.1e-11
    r = 14e-9
    x = d / 2

    def __init__(self, rozmiar_macierzy_blok):
        self.rozmiar_macierzy_blok = rozmiar_macierzy_blok
