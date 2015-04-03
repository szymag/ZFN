__author__ = 'szymag'

from src.Siec import Siec


class SiecTrojkatna(Siec):
    def __init__(self, dlugosc_a1, dlugosc_a2, kat, zakres_1, zakres_2):
        Siec.__init__(self)
        pass

    def wylicz_wspolczynniki(self, rodzaj):
        if rodzaj == 'okragly':
            return self.rdzen_okragly.wylicz_wspolczynniki(self.MoA, self.MoB, self.d, self.s, self.r)
        elif rodzaj == 'kwadratowy':
            return self.rdzen_kwadratowy.wylicz_wspolczynniki(self.MoA, self.MoB, self.d, self.s, self.r)