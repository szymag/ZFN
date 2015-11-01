__author__ = 'szymag'

import math

import numpy as np

from src.Rdzen import Rdzen


class RdzenKwadratowyTrojkatna(Rdzen):
    def __init__(self, lista_wektorow_b1, lista_wektorow_b2):
        Rdzen.__init__(self, lista_wektorow_b1, lista_wektorow_b2)

    def wzory_wspolczynniki_fouriera(self, wektor_g1, wektor_g2, MoA, MoB, d, s, r):
        wspolczynnik_fouriera = (MoA - MoB) * d ** 2 / s * \
                                np.sinc(np.array(wektor_g1 * d / 2 / math.pi)) * \
                                np.sinc(np.array(wektor_g2 * d / 2 / math.pi))
        return wspolczynnik_fouriera

    def wylicz_wspolczynniki_fouriera(self, MoA, MoB, d, s, r):
        wsp = self.tablica_wspolczynniki_fouriera(MoA, MoB, d, s, r)
        wsp[int(len(self.gx) / 2)][int(len(self.gy) / 2)][2] = \
            MoA * (d ** 2 / s) + MoB * (1 - d ** 2 / s)
        return wsp
