from numpy import pi


class ParametryMaterialowe:
    mu0H0 = 0.001
    H0 = mu0H0 / (4 * pi * 10 ** -7)
    d = 6
    a = 10
    s = a ** 2
    MoA = 10
    MoB = 4
    r = 3
    x = 1

    def __init__(self, rozmiar_macierzy_blok):
        self.rozmiar_macierzy_blok = rozmiar_macierzy_blok
