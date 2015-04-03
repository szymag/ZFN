__author__ = 'szymag'

import math

from src.Rdzen import Rdzen


class RdzenOkragly(Rdzen):
    def __init__(self, lista_wektorow_b1, lista_wektorow_b2):
        Rdzen.__init__(self, lista_wektorow_b1, lista_wektorow_b2)

    def wylicz_wspolczynniki(self, MoA, MoB, d, s, r, typ):
        """
        metoda zdefiniowana w klasach pochodnych
        :return: zwracana jest tablica zawierająca współczynniki fouriera dla rdzenia okrągłego.
        """
        wsp = self.wspolczynniki_fouriera_tablica(MoA, MoB, d, s, r, 'okragly')
        wsp[int(len(self.gx) / 2)][int(len(self.gy) / 2)][2] = \
            (MoA - MoB) * math.pi * r ** 2 / s + MoB
        return wsp