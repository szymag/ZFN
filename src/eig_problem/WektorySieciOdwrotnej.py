from itertools import product

from numpy import array, pi, sqrt


class WektorySieciOdwrotnej:
    def __init__(self, dlugosc_skladowej_y, dlugosc_skladowej_z, rozmiar_macierzy_blok):
        self.dlugosc_skladowej_y = dlugosc_skladowej_y
        self.dlugosc_skladowej_z = dlugosc_skladowej_z
        self.ilosc_wektorow = rozmiar_macierzy_blok

    def wspolrzedna_wektora(self, k):
        assert k == 1 or k == 2, 'k should bo 1 or 2'

        lista = array(range(int((-sqrt(self.ilosc_wektorow) + 1) / 2), int(sqrt(self.ilosc_wektorow) / 2) + 1))
        if k == 1:
            lista = [2 * pi * int(i) / self.dlugosc_skladowej_y for i in lista]
        elif k == 2:
            lista = [2 * pi * int(i) / self.dlugosc_skladowej_z for i in lista]

        return lista

    @property
    def lista_wektorow(self):
        lista = list(product(*(self.wspolrzedna_wektora(1), self.wspolrzedna_wektora(2))))
        return lista

    def wektory_odwrotne(self):
        lista1 = array(range(int(-self.ilosc_wektorow / 2), 0))
        lista1 = [((k * (k % 2)), (k * ((k + 1) % 2))) for k in lista1]
        lista2 = array(range(0, int(-self.ilosc_wektorow / 2)))
        lista2 = [((k * ((k + 1) % 2)), (k * (k % 2))) for k in lista2]
        return lista1 + lista2
