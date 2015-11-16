from itertools import product
from numpy import array, pi


class WektorySieciOdwrotnej:
    def __init__(self, dlugosc_skladowej_y, dlugosc_skladowej_z, ilosc_wektorow):
        self.dlugosc_skladowej_y = dlugosc_skladowej_y
        self.dlugosc_skladowej_z = dlugosc_skladowej_z
        self.ilosc_wektorow = ilosc_wektorow

    def wspolrzedna_wektora(self, k):
        assert k == 1 or k == 2, 'k should bo 1 or 2'

        lista = array(range(int(-self.ilosc_wektorow / 2), int(self.ilosc_wektorow / 2) + 1))
        if k == 1:
            lista = [2 * pi * int(i) / self.dlugosc_skladowej_y for i in lista]
        elif k == 2:
            lista = [2 * pi * int(i) / self.dlugosc_skladowej_z for i in lista]
        return lista

    @property
    def lista_wektorow(self):
        return list(product(*(self.wspolrzedna_wektora(1), self.wspolrzedna_wektora(2))))


