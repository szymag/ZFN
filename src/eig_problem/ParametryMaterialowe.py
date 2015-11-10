from numpy import pi


class ParametryMaterialowe:
    """
    :param d: grubość warstwy
    """
    mu0H0 = 0.001
    gamma = 194.6 * 10 ** 9
    H0 = mu0H0 / (4 * pi * 10 ** -7)
    d = 4 * 10 ** -9
    a = 30 * 10 ** -9
    MoA = 1.15 * 10 ** 6
    MoB = 0.558 * 10 ** 6
    r = 14 * 10 ** -9
    x = d / 2

    def __init__(self, rozmiar_macierzy_blok):
        self.rozmiar_macierzy_blok = rozmiar_macierzy_blok
