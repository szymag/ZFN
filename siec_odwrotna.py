__author__ = 'sz'
import math
import cmath

import numpy as np
import scipy.special


class SiatkaPunktow(object):
    def __init__(self, zakres_x, zakres_y, krok):
        self.zakres_x = zakres_x
        self.zakres_y = zakres_y
        self.krok = krok

    def generowanie_punktow(self, liczba):
        temp = 0
        lista = []
        while liczba >= temp:
            lista.extend([temp])
            temp += self.krok
        return lista

    def siatka(self):

        tempx = self.generowanie_punktow(self.zakres_x)
        tempy = self.generowanie_punktow(self.zakres_y)

        tablica = []
        for ii in tempx:
            temp = []
            for jj in tempy:
                temp.append([ii, jj, 0])
            tablica.append(temp)
        return tablica


class WektorySieci(object):
    a_z = [0, 0, 1]

    def __init__(self, dlugosc_a1, dlugosc_a2, kat, zakres_1, zakres_2):
        self.dlugosc_a1 = dlugosc_a1
        self.dlugosc_a2 = dlugosc_a2
        self.kat = kat
        self.zakres_1 = list(range(-zakres_1, zakres_1 + 1))
        self.zakres_2 = list(range(-zakres_2, zakres_2 + 1))

    def wektor_a1(self):
        return [self.dlugosc_a1, 0, 0]

    def wektor_a2(self):
        return [self.dlugosc_a2 * math.cos(math.radians(self.kat)),
                self.dlugosc_a2 * math.sin(math.radians(self.kat)),
                0]

    def wektor_b(self, k):
        x = self.wektor_a1()
        y = self.wektor_a2()
        z = self.a_z
        res = 2 * math.pi
        if k == 1:
            res *= np.cross(y, z)
        else:
            res *= np.cross(z, x)
        return res / np.dot(x, np.cross(y, z))

    def lista_wektorow_b1(self):
        return [self.wektor_b(1) * k for k in self.zakres_1]

    def lista_wektorow_b2(self):
        return [self.wektor_b(2) * k for k in self.zakres_2]

    def wektor_wypadkowy(self):
        return [a + b for a, b in zip(self.wektor_b(1), self.wektor_b(2))]


class SiecKwadratowa(WektorySieci):
    d = 6
    a = 10
    s = a ** 2
    MoA = 10
    MoB = 4
    R = 3

    def __init__(self, dlugosc_a1, dlugosc_a2, kat, zakres_1, zakres_2):
        WektorySieci.__init__(self, dlugosc_a1, dlugosc_a2, kat, zakres_1, zakres_2)

    def funkcja(self, k, liczba1, liczba2):

        if k == 1:
            wspolczynnik_fouriera = (self.MoA - self.MoB) * self.d ** 2 / self.s * \
                                    np.sinc(np.array(liczba1 * self.d / 2 / math.pi)) * \
                                    np.sinc(np.array(liczba2 * self.d / 2 / math.pi))
        else:
            wspolczynnik_fouriera = 2 * (self.MoA - self.MoB) * math.pi * self.R ** 2 / self.s * \
                                    scipy.special.j1(math.sqrt(liczba1 ** 2 + liczba2 ** 2) * self.R) \
                                    / (math.sqrt(liczba1 ** 2 + liczba2 ** 2) * self.R)
        return wspolczynnik_fouriera

    def rdzen_kwadratowy(self):
        gx = self.lista_wektorow_b1()
        gy = self.lista_wektorow_b2()
        wsp = []
        for ii in range(len(gx)):
            temp = []
            for jj in range(len(gy)):
                wspolczynnik_fouriera = self.funkcja(1, gx[ii][0], gy[jj][1])
                temp.append([gx[ii][0], gy[jj][1], wspolczynnik_fouriera])
            wsp.append(temp)
        wsp[int(len(gx) / 2)][int(len(gy) / 2)][2] = self.MoA * (self.d ** 2 / self.s) \
                                                     + self.MoB * (1 - self.d ** 2 / self.s)
        return wsp

    def rdzen_okragly(self):
        gx = self.lista_wektorow_b1()
        gy = self.lista_wektorow_b2()
        wsp = []
        for ii in range(len(gx)):
            temp = []
            for jj in range(len(gy)):
                wspolczynnik_fouriera = self.funkcja(2, gx[ii][0], gy[jj][1])
                temp.append([gx[ii][0], gy[jj][1], wspolczynnik_fouriera])
            wsp.append(temp)
        wsp[int(len(gx) / 2)][int(len(gy) / 2)][2] = \
            (self.MoA - self.MoB) * math.pi * self.R ** 2 / self.s + self.MoB
        return wsp


class Magnetyzacja(object):
    object1 = SiatkaPunktow(5, 5, 0.1)
    siatka = object1.siatka()
    object2 = SiecKwadratowa(10, 10, 90, 5, 5)
    kwadratowa_okragy = object2.rdzen_okragly()
    kwadratowa_kwadratowy = object2.rdzen_kwadratowy()

    def __init__(self):
        pass

    def magnetyzacja_w_punkcie(self, parametry, x, y):
        temp = 0
        for ii in range(len(parametry)):
            for jj in range(len(parametry)):
                temp += (parametry[ii][jj][2]
                         * (cmath.exp(complex(1j) * (parametry[ii][jj][0] * x + (parametry[ii][jj][1] * y))))).real
        return temp

    def magnetyzacja(self, parametry):
        table = self.siatka
        for ii in range(len(table)):
            for jj in range(len(table)):
                table[ii][jj][2] = \
                    self.magnetyzacja_w_punkcie(parametry, table[ii][jj][0], table[ii][jj][1])
        return table

    def magnetyzacja_dla_sieci(self, k):
        if k == 1:
            return self.magnetyzacja(self.kwadratowa_okragy)
        elif k == 2:
            return self.magnetyzacja(self.kwadratowa_kwadratowy)
        elif k == 3:
            return None
        elif k == 4:
            return None
        else:
            return None


w = Magnetyzacja()
print(w.magnetyzacja_dla_sieci(1))