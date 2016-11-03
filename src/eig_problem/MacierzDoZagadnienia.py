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

    def __init__(self, input_fft, ilosc_wektorow, wektor_q, a=ParametryMaterialowe.a, MoA=ParametryMaterialowe.MoA,
                 MoB=ParametryMaterialowe.MoB, lA=ParametryMaterialowe.lA,
                 lB=ParametryMaterialowe.lB, d=ParametryMaterialowe.d,
                 x=ParametryMaterialowe.x, H0=ParametryMaterialowe.H0):
        self.a = a
        self.ilosc_wektorow = ilosc_wektorow
        self.d = d
        self.x = x
        self.H0 = H0
        self.tmp = FFTfromFile(input_fft, ilosc_wektorow)
        self.macierz_M = np.zeros((2 * ilosc_wektorow, 2 * ilosc_wektorow), dtype=complex)
        self.magnetyzacja = self.tmp.fourier_coefficient(MoA, MoB)
        self.dlugosc_wymiany = self.tmp.fourier_coefficient(lA, lB)
        self.lista_wektorow = WektorySieciOdwrotnej(ilosc_wektorow).lista_wektorow2d('min')
        tmp = int(sqrt(self.ilosc_wektorow))
        self.shift = np.array([tmp - 1, tmp - 1])
        self.wektor_q = wektor_q

    def delta_kroneckera(self):
        """
        Metoda dodająca do odpowienich elementów macierzowych pierwszy element z wyrażenia na M: 1 lub -1
        """
        for i in range(self.ilosc_wektorow, 2 * self.ilosc_wektorow):
            self.macierz_M[i - self.ilosc_wektorow][i] += 1.
            self.macierz_M[i][i - self.ilosc_wektorow] -= 1.

    # noinspection PyTypeChecker
    def pole_wymiany_II(self, wektor_2):
        shift = self.shift
        a = self.a
        wektor_q = self.wektor_q
        H0 = self.H0
        tab_from_wektor_1 = np.asfortranarray(
            np.broadcast_to(self.lista_wektorow, (self.ilosc_wektorow, self.ilosc_wektorow, 2)), dtype=int)
        tab_from_vec_l = np.transpose(tab_from_wektor_1, (1, 0, 2))
        tmp = ne.evaluate('tab_from_wektor_1 - tab_from_vec_l + shift')
        tmp1 = self.magnetyzacja[tmp[:, :, 0], tmp[:, :, 1]]
        tmp2 = ne.evaluate('tab_from_vec_l - wektor_2 + shift')
        tmp3 = self.dlugosc_wymiany[tmp2[:, :, 0], tmp2[:, :, 1]]
        tmp4 = np.dot(ne.evaluate('2 * 3.14159265 * tab_from_vec_l / a + wektor_q'),
                      ne.evaluate('2 * 3.14159265 * wektor_2 / a + wektor_q'))
        return ne.evaluate('sum(tmp1*tmp3*tmp4/H0, 0)')

    # noinspection PyTypeChecker
    def trzecie_wyrazenie_xy(self, wektor_2):
        """
        Metoda obliczająca człon elmentu macierzowego pochodzenia dipolowego.
        :param wektor_2: j-ty wektor.
        :return: Wynikiem jest drugi wyraz sumy.
        """
        H0 = self.H0
        wekt_2 = np.array(2 * np.pi * wektor_2 / self.a)
        norm = sqrt((self.wektor_q[0] + wekt_2[0]) ** 2 + (self.wektor_q[1] + wekt_2[1]) ** 2)
        tmp1 = (self.wektor_q[1] + wekt_2[1]) ** 2 / norm ** 2
        tmp2 = 1 - cosh(norm * self.x) * exp(-norm * self.d / 2.)
        tmp3 = self.lista_wektorow - wektor_2 + self.shift
        tmp4 = self.magnetyzacja[tmp3[:, 0], tmp3[:, 1]]
        return ne.evaluate('tmp1 * tmp2 * tmp4 / H0')

    def trzecie_wyrazenie_yx(self, wektor_2):
        H0 = self.H0
        wekt_2 = np.array(2 * np.pi * wektor_2 / self.a)
        norm = sqrt((self.wektor_q[0] + wekt_2[0]) ** 2 + (self.wektor_q[1] + wekt_2[1]) ** 2)
        tmp1 = cosh(norm * self.x) * exp(-norm * self.d / 2.)
        tmp3 = self.lista_wektorow - wektor_2 + self.shift
        tmp4 = self.magnetyzacja[tmp3[:, 0], tmp3[:, 1]]
        return ne.evaluate('tmp1 * tmp4 / H0')

    # noinspection PyTypeChecker
    def czwarte_wyrazenie(self, wektor_2):
        """
        Metoda obliczająca człon elmentu macierzowego pochodzenia dipolowego.
        :param wektor_2: j-ty wektor.
        :return: Wynikiem jest czwarte wyrażenie w sumie na element macierzy M.
        """
        H0 = self.H0
        wektor_1 = self.lista_wektorow
        wekt_2 = np.array(2 * np.pi * wektor_2 / self.a)
        wekt_1 = 2 * np.pi * wektor_1 / self.a
        tmp1 = (wekt_1[:, 0] - wekt_2[0]) ** 2 / (np.linalg.norm(wekt_1 - wekt_2, axis=1) ** 2 + 1e-36)
        tmp2 = 1 - np.cosh(np.linalg.norm(wekt_1 - wekt_2, axis=1) * self.x) * \
                np.exp(-np.linalg.norm(wekt_1 - wekt_2, axis=1) * self.d / 2.)
        tmp3 = wektor_1 - wektor_2 + self.shift
        tmp4 = self.magnetyzacja[tmp3[:, 0], tmp3[:, 1]]
        return ne.evaluate('tmp1 * tmp2 * tmp4 / H0')

    def wypelnienie_macierzy(self):
        """
        Główna metoda tej klasy. Wywołuje ona dwie metody: 'macierz_xy' oraz 'macierz_yx. W pętli, dla każdego elementu
        z odpowiednich macierzy blokowych wypełnia je.
        :return: Tablica do zagadnienia własnego.
        """
        self.delta_kroneckera()
        w2 = self.lista_wektorow
        pool = Pool()
        tmp1 = np.array(pool.map(self.pole_wymiany_II, w2))
        #tmp1 = np.apply_along_axis(self.pole_wymiany_II, 1, w2)
        tmp4 = np.array(pool.map(self.czwarte_wyrazenie, w2))
        tmp2 = np.array(pool.map(self.trzecie_wyrazenie_xy, w2))
        tmp3 = np.array(pool.map(self.trzecie_wyrazenie_yx, w2))
        self.macierz_M[self.ilosc_wektorow:, 0:self.ilosc_wektorow] += -tmp1 - tmp3 + tmp4
        self.macierz_M[0:self.ilosc_wektorow, self.ilosc_wektorow:] += tmp1 + tmp2 - tmp4
        pool.close()
        return self.macierz_M

    def wypisz_macierz(self):
        """
        :return: Wypisuje tablice do pliku tekstowego.
         Ważne! Przed wypisaniem, należy wypełnić macierz_M metodą 'wypełnienie_macierzy'
        """
        self.wypelnienie_macierzy()
        np.savetxt('macierz.txt', np.array(self.macierz_M))


if __name__ == "__main__":
    q = MacierzDoZagadnienia('radius90.txt', 9, np.array([1e-9, 0]))
    q.wypisz_macierz()
