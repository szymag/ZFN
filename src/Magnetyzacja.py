__author__ = 'szymag'

import cmath
from itertools import product

class Magnetyzacja(object):
    def __init__(self, siatka_punktow, siec):
        self.siatka_punktow = siatka_punktow
        self.siec = siec

    def magnetyzacja_w_punkcie(self, parametry, x, y):
        temp = 0
        for ii, jj in product(range(len(parametry)), range(len(parametry))):
                temp += (parametry[ii][jj][2]
                         * (cmath.exp(complex(1j) * (parametry[ii][jj][0] * x + (parametry[ii][jj][1] * y))))).real
        return temp

    def magnetyzacja(self, parametry):
        table = self.siatka_punktow.siatka()
        for ii, jj in product(range(len(table)), range(len(table[0]))):
            table[ii][jj][2] = self.magnetyzacja_w_punkcie(parametry, table[ii][jj][0], table[ii][jj][1])
        return table

    def magnetyzacja_dla_sieci(self, typ_rdzenia):
        if typ_rdzenia == 'kwadratowy':
            return self.magnetyzacja(self.siec.wylicz_wspolczynniki_fouriera('kwadratowy'))
        elif typ_rdzenia == 'okragly':
            return self.magnetyzacja(self.siec.wylicz_wspolczynniki_fouriera('okragly'))
        else:
            return None

    def magnetyzacja1(self, parametry):
        lista = self.siatka_punktow.siatka()[0]
        for ii in range(len(lista)):
            lista[ii][2] = self.magnetyzacja_w_punkcie(parametry, lista[ii][0], lista[ii][1])
        return lista

    def magnetyzacja_pod_plot(self, typ_rdzenia):
        lista = self.siec
        if typ_rdzenia == 'kwadratowy':
            return self.magnetyzacja1(lista.wylicz_wspolczynniki_fouriera('kwadratowy'))
        elif typ_rdzenia == 'okragly':
            return self.magnetyzacja1(lista.wylicz_wspolczynniki_fouriera('okragly'))
        else:
            return None

