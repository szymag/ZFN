__author__ = 'szymag'

import math

import numpy as np
import scipy.special


class Rdzen:
    def __init__(self, lista_wektorow_b1, lista_wektorow_b2):
        self.gx = lista_wektorow_b1[:]
        self.gy = lista_wektorow_b2[:]

    def wspolczynnik_fouriera(self, k, liczba1, liczba2, MoA, MoB, d, s, R):
        if k == 1:
            wspolczynnik_fouriera = (MoA - MoB) * d ** 2 / s * \
                                    np.sinc(np.array(liczba1 * d / 2 / math.pi)) * \
                                    np.sinc(np.array(liczba2 * d / 2 / math.pi))
        else:
            wspolczynnik_fouriera = 2 * (MoA - MoB) * math.pi * R ** 2 / s * \
                                    scipy.special.j1(math.sqrt(liczba1 ** 2 + liczba2 ** 2) * R) \
                                    / (math.sqrt(liczba1 ** 2 + liczba2 ** 2) * R)
        return wspolczynnik_fouriera

    def wylicz_wspolczynniki_fouriera(self, MoA, MoB, d, s, R, k):
        wsp = []
        for ii in range(len(self.gx)):
            temp = []
            for jj in range(len(self.gy)):
                wspolczynnik_fouriera = self.wspolczynnik_fouriera(k, self.gx[ii][0], self.gy[jj][1], MoA, MoB, d, s, R)
                temp.append([self.gx[ii][0], self.gy[jj][1], wspolczynnik_fouriera])
            wsp.append(temp)
        return wsp

    def wylicz_wspolczynniki(self, MoA, MoB, d, s, R, k):
        """ metoda zdefiniowana w klasach pochodnych """
        pass