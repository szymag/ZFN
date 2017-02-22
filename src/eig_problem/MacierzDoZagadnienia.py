# -*- coding: utf-8 -*-
from math import sqrt, exp, cosh
import numpy as np
from src.eig_problem.FFTfromFile import FFTfromFile
from src.eig_problem.ParametryMaterialowe import ParametryMaterialowe
from src.eig_problem.WektorySieciOdwrotnej import WektorySieciOdwrotnej
from multiprocessing import Pool
import numexpr as ne


class MacierzDoZagadnienia:
    """
    Klasa, w której tworzona jest macierz zagadnienia własnego. Pobiera ona współczynniki Fouriera z klasy FFTfromFile.
    Współczynniki są dla niej tworzone poprzez FFT z pliku graficznego.
    """

    def __init__(self, input_fft, rec_vector_x, rec_vector_y, wektor_q, a=ParametryMaterialowe.a, b=ParametryMaterialowe.b,
                 MoA=ParametryMaterialowe.MoA, MoB=ParametryMaterialowe.MoB, lA=ParametryMaterialowe.lA,
                 lB=ParametryMaterialowe.lB, d=ParametryMaterialowe.d,
                 x=ParametryMaterialowe.x, H0=ParametryMaterialowe.H0):
        self.a = a
        self.b = b
        self.d = d
        self.rec_vector_x = rec_vector_x
        self.rec_vector_y = rec_vector_y
        self.ilosc_wektorow = rec_vector_x * rec_vector_y
        self.x = x
        self.H0 = H0
        self.tmp = FFTfromFile(input_fft, (2*self.rec_vector_x-1, 2*self.rec_vector_y-1))
        self.macierz_M = np.zeros((2 * self.ilosc_wektorow, 2 * self.ilosc_wektorow), dtype=complex)
        self.magnetyzacja = self.tmp.fourier_coefficient(MoA, MoB)
        self.dlugosc_wymiany = self.tmp.fourier_coefficient(lA, lB)
        self.lista_wektorow = WektorySieciOdwrotnej(self.ilosc_wektorow).lista_wektorow2d('min')
        self.shift = np.array([self.rec_vector_x - 1, self.rec_vector_y - 1])
        self.wektor_q = wektor_q

    def save_to_file_matrix(self):
        self.fill_matrix()
        np.savetxt('macierz.txt', np.array(self.macierz_M))

    def fill_matrix(self):
        self.kroneker_delta()
        w2 = self.lista_wektorow
        pool = Pool()
        tmp1 = np.array(pool.map(self.exchange_field, w2))
        tmp4 = np.array(pool.map(self.static_demagnetizing_field, w2))
        tmp2 = np.array(pool.map(self.dynamic_demagnetizing_field_in_plane, w2))
        tmp3 = np.array(pool.map(self.dynamic_demagnetizing_field_out_of_plane, w2))
        self.macierz_M[self.ilosc_wektorow:, 0:self.ilosc_wektorow] += -tmp1 - tmp3 + tmp4
        self.macierz_M[0:self.ilosc_wektorow, self.ilosc_wektorow:] += tmp1 + tmp2 - tmp4
        pool.close()
        return self.macierz_M

    def kroneker_delta(self):
        for i in range(self.rec_vector_x, 2 * self.rec_vector_x):
            self.macierz_M[i - self.rec_vector_x][i] += 1.
            self.macierz_M[i][i - self.rec_vector_x] -= 1.

    # noinspection PyTypeChecker
    def exchange_field(self, wektor_2):
        shift = self.shift
        a = self.a
        b = self.b
        wektor_q = self.wektor_q
        H0 = self.H0
        tab_from_wektor_1 = np.asfortranarray(
            np.broadcast_to(self.lista_wektorow, (self.ilosc_wektorow , self.ilosc_wektorow, 2)), dtype=int)
        tab_from_vec_l = np.transpose(tab_from_wektor_1, (1, 0, 2))
        tmp = ne.evaluate('tab_from_wektor_1 - tab_from_vec_l + shift')
        tmp1 = self.magnetyzacja[tmp[:, :, 0], tmp[:, :, 1]]
        tmp2 = ne.evaluate('tab_from_vec_l - wektor_2 + shift')
        tmp3 = self.dlugosc_wymiany[tmp2[:, :, 0], tmp2[:, :, 1]]
        tab_from_vec_l = tab_from_vec_l / [b, a]
        wektor_2 = wektor_2 / [b, a]
        tmp4 = np.dot(ne.evaluate('2 * 3.14159265 * tab_from_vec_l + wektor_q'),
                      ne.evaluate('2 * 3.14159265 * wektor_2 + wektor_q'))
        return ne.evaluate('sum(tmp1*tmp3*tmp4/H0, 0)')

    # noinspection PyTypeChecker
    def dynamic_demagnetizing_field_in_plane(self, wektor_2):
        H0 = self.H0
        wekt_2 = np.array(2 * np.pi * wektor_2 / [self.b, self.a])
        norm = sqrt((self.wektor_q[0] + wekt_2[0]) ** 2 + (self.wektor_q[1] + wekt_2[1]) ** 2)
        tmp1 = (self.wektor_q[0] + wekt_2[0]) ** 2 / norm ** 2
        tmp2 = 1 - cosh(norm * self.x) * exp(-norm * self.d / 2.)
        tmp3 = self.lista_wektorow - wektor_2 + self.shift
        tmp4 = self.magnetyzacja[tmp3[:, 0], tmp3[:, 1]]
        return ne.evaluate('tmp1 * tmp2 * tmp4 / H0')

    # noinspection PyTypeChecker
    def dynamic_demagnetizing_field_out_of_plane(self, wektor_2):
        H0 = self.H0
        wekt_2 = np.array(2 * np.pi * wektor_2 / [self.b, self.a])
        norm = sqrt((self.wektor_q[0] + wekt_2[0]) ** 2 + (self.wektor_q[1] + wekt_2[1]) ** 2)
        tmp1 = cosh(norm * self.x) * exp(-norm * self.d / 2.)
        #print(tmp1)
        tmp3 = self.lista_wektorow - wektor_2 + self.shift
        tmp4 = self.magnetyzacja[tmp3[:, 0], tmp3[:, 1]]
        #print(tmp3)
        return ne.evaluate('tmp1 * tmp4 / H0')

    # noinspection PyTypeChecker
    def static_demagnetizing_field(self, wektor_2):
        # along external field
        H0 = self.H0
        wektor_1 = self.lista_wektorow
        wekt_2 = np.array(2 * np.pi * wektor_2 / [self.b, self.a])
        wekt_1 = 2 * np.pi * wektor_1 / [self.b, self.a]
        tmp1 = (wekt_1[:, 1] - wekt_2[1]) ** 2 / (np.linalg.norm(wekt_1 - wekt_2, axis=1) ** 2 + 1e-36)
        tmp2 = 1 - np.cosh(np.linalg.norm(wekt_1 - wekt_2, axis=1) * self.x) * \
                   np.exp(-np.linalg.norm(wekt_1 - wekt_2, axis=1) * self.d / 2.)
        tmp3 = wektor_1 - wektor_2 + self.shift
        tmp4 = self.magnetyzacja[tmp3[:, 0], tmp3[:, 1]]
        #return ne.evaluate('tmp1 * tmp2 * tmp4 / H0')
        return tmp4 * tmp1 * tmp2 /self.H0

if __name__ == "__main__":
    q = MacierzDoZagadnienia('ff=0.5.txt', 3, 3, np.array([1e-9, 0]))
    q.save_to_file_matrix()
