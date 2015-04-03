__author__ = 'szymag'

import math

import numpy as np
import scipy.special


class Rdzen:
    """
    Klasa. w której obliczane są współczynniki fourierea, dla różnych typów sieci.
    """
    def __init__(self, lista_wektorow_b1, lista_wektorow_b2):
        """
        :param lista_wektorow_b1: definiowanie listy wektorów sieci odwrtonej wzdłuż osi OX
        :param lista_wektorow_b2: definiowanie kolejnej listy wektorów sieci odwrotnej
        """
        self.gx = lista_wektorow_b1[:]
        self.gy = lista_wektorow_b2[:]


    def wspolczynniki_fouriera(self, typ_rdzenia, wektor_g1, wektor_g2, MoA, MoB, d, s, r):
        """
        :param typ_rdzenia: oznacza, jaki typ rdzenia zostaje wprowadzony.
        dla typ = 'kwadratowy': rdzeń jest kwadratowy
        dla typ = 'okragly': rdzeń jest okrągły
        :param wektor_g1:
        :param wektor_g2:
        :param MoA: magnetyzacja rdzenia
        :param MoB: magnetyzacja wypełnienia
        :param d: bok rdzenia kwadratowego
        :param s: powirzchnia komórki elementarnej
        :param r: promień rdzenia
        :return: zwracany jest współczynnik fouriera
        """
        global wspolczynnik_fouriera
        if typ_rdzenia == 'kwadratowy':
            wspolczynnik_fouriera = (MoA - MoB) * d ** 2 / s * \
                                    np.sinc(np.array(wektor_g1 * d / 2 / math.pi)) * \
                                    np.sinc(np.array(wektor_g2 * d / 2 / math.pi))
        elif typ_rdzenia == 'okragly':
            wspolczynnik_fouriera = 2 * (MoA - MoB) * math.pi * r ** 2 / s * \
                                    scipy.special.j1(math.sqrt(wektor_g1 ** 2 + wektor_g2 ** 2) * r) \
                                    / (math.sqrt(wektor_g1 ** 2 + wektor_g2 ** 2 + (10 ** -10)) * r)
        return wspolczynnik_fouriera

    def wspolczynniki_fouriera_tablica(self, MoA, MoB, d, s, r, typ_rdzenia):
        """
        :param MoA: magnetyzacja rdzenia
        :param MoB: magnetyzacja wypełnienia
        :param d: bok rdzenia kwadratowego
        :param s: powirzchnia komórki elementarnej
        :param r: promień rdzenia
        :return: zwracana jest tablica zawierająca wszystkie współczynniki fourieria dla dla par wektor_g1, wektor_g2
        """
        wsp = []
        for ii in range(len(self.gx)):
            temp = []
            for jj in range(len(self.gy)):
                wspolczynnik_fouriera = self.wspolczynniki_fouriera(typ_rdzenia, self.gx[ii][0], self.gy[jj][1],
                                                                    MoA, MoB, d, s, r)
                temp.append([self.gx[ii][0], self.gy[jj][1], wspolczynnik_fouriera])
            wsp.append(temp)
        return wsp

    def wylicz_wspolczynniki(self, MoA, MoB, d, s, r, typ_sieci):
        """ metoda zdefiniowana w klasach pochodnych """
        pass