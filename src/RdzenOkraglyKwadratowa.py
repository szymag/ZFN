__author__ = 'szymag'

import math

import scipy.special

from src.Rdzen import Rdzen


class RdzenOkraglyKwadratowa(Rdzen):
    def __init__(self, lista_wektorow_b1, lista_wektorow_b2):
        Rdzen.__init__(self, lista_wektorow_b1, lista_wektorow_b2)

    def wzory_wspolczynniki_fouriera(self, wektor_g1, wektor_g2, MoA, MoB, d, s, r):
        wspolczynnik_fouriera = 2 * (MoA - MoB) * math.pi * r ** 2 / s * \
                                scipy.special.j1(math.sqrt(wektor_g1 ** 2 + wektor_g2 ** 2) * r) \
                                / (math.sqrt(wektor_g1 ** 2 + wektor_g2 ** 2 + (10 ** -10)) * r)
        return wspolczynnik_fouriera

    def wylicz_wspolczynniki_fouriera(self, MoA, MoB, d, s, r):
        wsp = self.tablica_wspolczynniki_fouriera(MoA, MoB, d, s, r)
        wsp[int(len(self.gx) / 2)][int(len(self.gy) / 2)][2] = \
            (MoA - MoB) * math.pi * r ** 2 / s + MoB
        return wsp