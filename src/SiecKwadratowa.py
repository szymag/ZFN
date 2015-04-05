__author__ = 'szymag'

from src.Siec import Siec
from src.WektorySieci import WektorySieci
from src.RdzenOkragly import RdzenOkragly
from src.RdzenKwadratowy import RdzenKwadratowy


class SiecKwadratowa(Siec):
    def __init__(self, dlugosc_a1, dlugosc_a2, zakres_wektorow_gx, zakres_wektorow_gy):
        Siec.__init__(self, typ_sieci='kwadratowa')
        self.wektory_sieci = WektorySieci(dlugosc_a1, dlugosc_a2, 90, zakres_wektorow_gx, zakres_wektorow_gy)
        self.rdzen_okragly = RdzenOkragly(self.wektory_sieci.lista_wektorow_b1(),
            self.wektory_sieci.lista_wektorow_b2(), typ_sieci='kwadratowa')
        self.rdzen_kwadratowy = RdzenKwadratowy(self.wektory_sieci.lista_wektorow_b1(),
            self.wektory_sieci.lista_wektorow_b2(), typ_sieci='kwadratowa')


    def wylicz_wspolczynniki(self, typ_rdzenia):
        if typ_rdzenia == 'okragly':
            return self.rdzen_okragly.wylicz_wspolczynniki(self.MoA, self.MoB, self.d, self.s, self.r, self.typ_sieci)
        elif typ_rdzenia == 'kwadratowy':
            return self.rdzen_kwadratowy.wylicz_wspolczynniki(self.MoA, self.MoB, self.d, self.s, self.r, 'kwadratowa')
