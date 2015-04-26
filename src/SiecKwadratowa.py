__author__ = 'szymag'

from src.Siec import Siec
from src.WektorySieci import WektorySieci
from src.RdzenOkraglyKwadratowa import RdzenOkraglyKwadratowa
from src.RdzenKwadratowyKwadratowa import RdzenKwadratowyKwadratowa


class SiecKwadratowa(Siec):
    def __init__(self, dlugosc_a1, dlugosc_a2, zakres_wektorow_gx, zakres_wektorow_gy):
        Siec.__init__(self, typ_sieci='kwadratowa')
        self.wektory_sieci = WektorySieci(dlugosc_a1, dlugosc_a2, 90, zakres_wektorow_gx, zakres_wektorow_gy)
        self.rdzen_okragly = RdzenOkraglyKwadratowa(self.wektory_sieci.lista_wektorow_b1(),
            self.wektory_sieci.lista_wektorow_b2())
        self.rdzen_kwadratowy = RdzenKwadratowyKwadratowa(self.wektory_sieci.lista_wektorow_b1(),
            self.wektory_sieci.lista_wektorow_b2())

    def wylicz_wspolczynniki_fouriera(self, typ_rdzenia):
        if typ_rdzenia == 'okragly':
            return self.rdzen_okragly.wylicz_wspolczynniki_fouriera(self.MoA, self.MoB, self.d, self.s, self.r)
        elif typ_rdzenia == 'kwadratowy':
            return self.rdzen_kwadratowy.wylicz_wspolczynniki_fouriera(self.MoA, self.MoB, self.d, self.s, self.r)
