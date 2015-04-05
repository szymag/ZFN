__author__ = 'szymag'

from src.Rdzen import Rdzen


class RdzenOkraglyKwadratowa(Rdzen):
    def __init__(self, lista_wektorow_b1, lista_wektorow_b2, typ_sieci):
        Rdzen.__init__(self, lista_wektorow_b1, lista_wektorow_b2, typ_sieci)

    def wylicz_wspolczynniki(self, MoA, MoB, d, s, r, typ_sieci):
        Rdzen.wspolczynniki_fouriera_tablica(self, )