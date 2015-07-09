__author__ = 'szymag'

from src.Siec import Siec
from src.WektorySieci import WektorySieci
from src.RdzenOkraglyTrojkatna import RdzenOkraglyTrojkatna
from src.RdzenKwadratowyTrojkatna import RdzenKwadratowyTrojkatna


class SiecTrojkatna(Siec):
    def __init__(self, dlugosc_a1, dlugosc_a2, zakres_wektorow_gx, zakres_wektorow_gy):
        Siec.__init__(self, typ_sieci='trojkatna')
        self.wektory_sieci = WektorySieci(dlugosc_a1, dlugosc_a2, 60, zakres_wektorow_gx, zakres_wektorow_gy)
        self.rdzen_okragly = RdzenOkraglyTrojkatna(self.wektory_sieci.lista_wektorow_b1(),
            self.wektory_sieci.lista_wektorow_b2())
        self.rdzen_kwadratowy = RdzenKwadratowyTrojkatna(self.wektory_sieci.lista_wektorow_b1(),
            self.wektory_sieci.lista_wektorow_b2())

    def wylicz_wspolczynniki_fouriera(self, typ_rdzenia):
        if typ_rdzenia == 'okragly':
            return self.rdzen_okragly.wylicz_wspolczynniki_fouriera(self.MoA, self.MoB, self.d, self.s, self.r)
        elif typ_rdzenia == 'kwadratowy':
            return self.rdzen_kwadratowy.wylicz_wspolczynniki_fouriera(self.MoA, self.MoB, self.d, self.s, self.r)
