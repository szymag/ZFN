from math import pi, sqrt

import numpy as np


class WektorySieciOdwrotnej:

    def __init__(self, rozmiar_macierzy_blok):
        self.ilosc_wektorow = rozmiar_macierzy_blok

    def lista_wektorow2d(self, typ):

        assert sqrt(self.ilosc_wektorow) == int(sqrt(self.ilosc_wektorow)) and self.ilosc_wektorow % 2 != 0, \
            'size of bolck matrix shuld gave natural number of sqrt and be odd'
        assert typ == 'max' or typ == 'min', 'typ should by max or min'
        if typ == 'max':
            indeks = int(sqrt(self.ilosc_wektorow)) - 1
            a = np.arange(-indeks, indeks +1)
            return np.array(np.meshgrid(a, a)).T.reshape(-1,2)
        else:
            indeks = int(sqrt(self.ilosc_wektorow + 1) / 2)
            a = np.arange(-indeks, indeks + 1)
            return np.array(np.meshgrid(a, a)).T.reshape(-1,2)

    def wspolrzedna_wektora1d(self, typ):
        assert typ == 'max' or typ == 'min', 'typ should by max or min'
        if typ == 'min':
            return np.arange(-(self.ilosc_wektorow - 1) // 2, (self.ilosc_wektorow - 1) // 2 + 1)
        elif typ == 'max':
            return np.arange(- (self.ilosc_wektorow - 1), (self.ilosc_wektorow - 1) + 1)

    def lista_wektorow1d(self, typ):
        return self.wspolrzedna_wektora1d(typ)

if __name__ == "__main__":
    q = WektorySieciOdwrotnej(9).lista_wektorow2d('max')
    print(q/[3,4])
