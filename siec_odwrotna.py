__author__ = 'sz'
import math

import numpy as np
import scipy.special


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


# q = WektorySieci(10, 10, 90, 5, 5)

# a = q.lista_wektorow_b1()
# print(a[0][0])


class SiecKwadratowa(WektorySieci):
    d = 6
    a = 10
    s = a ** 2
    MoA = 10
    MoB = 4
    R = 3

    def __init__(self, dlugosc_a1, dlugosc_a2, kat, zakres_1, zakres_2):
        WektorySieci.__init__(self, dlugosc_a1, dlugosc_a2, kat, zakres_1, zakres_2)
        self.gx = self.lista_wektorow_b1()
        self.gy = self.lista_wektorow_b2()


    def rdzen_kwadratowy(self):
        s = self.s
        d = self.d
        wsp = []
        for ii in range(len(self.gx)):
            temp = []
            for jj in range(len(self.gy)):
                wspolczynnik_fouriera = (self.MoA - self.MoB) * d ** 2 / s * \
                                        np.sinc(np.array(self.gx[ii][0] * d / 2 / math.pi)) * np.sinc(
                    np.array(self.gy[jj][1] * d / 2 / math.pi))
                temp.append([self.gx[ii][0], self.gy[jj][1], wspolczynnik_fouriera])
            wsp.append(temp)
        wsp[int(len(self.gx) / 2)][int(len(self.gy) / 2)][2] = self.MoA * (d ** 2 / s) + self.MoB * (1 - d ** 2 / s)
        return wsp

    def rdzen_okragly(self, funkcja):
        s = self.s
        gx = self.lista_wektorow_b1()
        gy = self.lista_wektorow_b2()
        MoA = self.MoA
        MoB = self.MoB
        R = self.R
        wsp = []
        # wsp = [ [funkcja(jj) for jj in range(len(gy))] for i in range(len(gx))]
        for ii in range(len(gx)):
            temp = []

            for jj in range(len(gy)):
                wspolczynnik_fouriera = 2 * (MoA - MoB) * math.pi * R ** 2 / s * \
                                        scipy.special.j1(math.sqrt(gx[ii][0] ** 2 + gy[jj][1] ** 2) * R) / (
                                            math.sqrt(gx[ii][0] ** 2 + gy[jj][1] ** 2) * R)
                temp.append([gx[ii][0], gy[jj][1], wspolczynnik_fouriera])
            wsp.append(temp)
        wsp[int(len(gx) / 2)][int(len(gy) / 2)][2] = (MoA - MoB) * math.pi * R ** 2 / s + MoB
        return wsp


q = SiecKwadratowa(10, 10, 90, 5, 5)
print(q.rdzen_okragly())
print("")
print(q.rdzen_kwadratowy())