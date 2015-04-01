__author__ = 'szymag'

from src.Rdzen import Rdzen


class RdzenKwadratowy(Rdzen):
    def __init__(self, lista_wektorow_b1, lista_wektorow_b2):
        Rdzen.__init__(self, lista_wektorow_b1, lista_wektorow_b2)

    def wylicz_wspolczynniki(self, MoA, MoB, d, s, R):
        """ metoda zdefiniowana w klasach pochodnych """
        wsp = self.wylicz_wspolczynniki_fouriera(MoA, MoB, d, s, R)
        wsp[int(len(self.gx) / 2)][int(len(self.gy) / 2)][2] = \
            MoA * (d ** 2 / s) + MoB * (1 - d ** 2 / s)
        return wsp

