__author__ = 'szymag'

from src.Siec import Siec
from src.WektorySieci import WektorySieci
from src.RdzenOkragly import RdzenOkragly
from src.RdzenKwadratowy import RdzenKwadratowy


class SiecKwadratowa(Siec):
    def __init__(self, dlugosc_a1, dlugosc_a2, kat, zakres_wektorow_gx, zakres_wektorow_gy):
        Siec.__init__(self)
        self.wektory_sieci = WektorySieci(dlugosc_a1, dlugosc_a2, kat, zakres_wektorow_gx, zakres_wektorow_gy)
        self.rdzen_okragly = RdzenOkragly(self.wektory_sieci.lista_wektorow_b1(),
            self.wektory_sieci.lista_wektorow_b2())
        self.rdzen_kwadratowy = RdzenKwadratowy(self.wektory_sieci.lista_wektorow_b1(),
            self.wektory_sieci.lista_wektorow_b2())

    def wylicz_wspolczynniki(self, rodzaj):
        if rodzaj == 'okragly':
            return self.rdzen_okragly.wylicz_wspolczynniki(self.MoA, self.MoB, self.d, self.s, self.R, 'okragly')
        elif rodzaj == 'kwadratowy':
            return self.rdzen_kwadratowy.wylicz_wspolczynniki(self.MoA, self.MoB, self.d, self.s, self.R, 'kwadratowy')
