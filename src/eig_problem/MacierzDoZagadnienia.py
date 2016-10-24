# -*- coding: utf-8 -*-
from src.eig_problem.cProfiler import do_cprofile
from math import exp, cosh, sqrt
import numpy as np
from src.eig_problem.FFTfromFile import FFTfromFile
from src.eig_problem.ParametryMaterialowe import ParametryMaterialowe
from src.eig_problem.WektorySieciOdwrotnej import WektorySieciOdwrotnej
from multiprocessing import Pool


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
        tmp = int(sqrt(len(self.lista_wektorow)))
        self.shift = np.array([tmp - 1, tmp - 1])
        self.wektor_q = wektor_q

    @staticmethod
    def norma_wektorow(wektor_1, wektor_2, znak):
        if znak == "+":
            return np.linalg.norm(wektor_1 + wektor_2)
        else:
            return np.linalg.norm(wektor_1 - wektor_2)

    def funkcja_c(self, wektor_1, wektor_2, znak):
        """
        Metoda obliczająca wartość funkcji C zdefinoweanej wzorem: f(g, x) = cosh(|g|x)*exp(-|g|d/2), gdzie g jest
        wektorem, a x współrzędną iksową, tzn. miejscem na warstwie, dla którego wykreśla się zależność dyspersyjną.
        :type znak: str
        :param wektor_1: Pierwszy wektor do obliczenia różnicy.
        :param wektor_2: Drugi wektor do obliczenia różnicy.
        :param znak: określa czy obliczana ma być różnica, czy suma wektorów.
        :return: wartość funkcji C.
        """
        wekt_wypadkowy = self.norma_wektorow(wektor_1, wektor_2, znak)
        return cosh(wekt_wypadkowy * self.x) * exp(-wekt_wypadkowy * self.d / 2.)

    def delta_kroneckera(self):
        """
        Metoda dodająca do odpowienich elementów macierzowych pierwszy element z wyrażenia na M: 1 lub -1
        """
        for i in range(self.ilosc_wektorow, 2 * self.ilosc_wektorow):
            self.macierz_M[i - self.ilosc_wektorow][i] += 1.
            self.macierz_M[i][i - self.ilosc_wektorow] -= 1.

    # noinspection PyTypeChecker
    def pole_wymiany_II(self, wektor_2):
        tab_from_wektor_1 = np.broadcast_to(self.lista_wektorow, (self.ilosc_wektorow, self.ilosc_wektorow, 2))
        tab_from_vec_l = np.transpose(tab_from_wektor_1, (1, 0, 2))
        tmp = tab_from_vec_l - tab_from_wektor_1 + self.shift
        tmp1 = self.magnetyzacja[tmp[:, :, 0], tmp[:, :, 1]]
        tmp2 = tab_from_vec_l - wektor_2 + self.shift
        tmp3 = self.dlugosc_wymiany[tmp2[:, :, 0], tmp2[:, :, 1]]
        tmp4 = np.dot(2 * np.pi * tab_from_vec_l / self.a + self.wektor_q,
                      2 * np.pi * wektor_2 / self.a + self.wektor_q)
        return np.sum(tmp1 * tmp3 * tmp4, axis=0) / self.H0

    # noinspection PyTypeChecker
    def trzecie_wyrazenie_xy(self, wektor_2):
        """
        Metoda obliczająca człon elmentu macierzowego pochodzenia dipolowego.
        :param wektor_1: i-ty wektor.
        :param wektor_2: j-ty wektor.
        :param wektor_q: Blochowski wektor. Jest on "uciąglony". Jest on zmienną przy wyzn aczaniu dyspersji.
        :return: Wynikiem jest drugi wyraz sumy.
        """
        wekt_2 = 2 * np.pi * wektor_2 / self.a
        tmp1 = (self.wektor_q[1] + wekt_2[1]) ** 2 / \
               ((self.wektor_q[0] + wekt_2[0]) ** 2 + (self.wektor_q[1] + wekt_2[1]) ** 2)
        tmp2 = 1 - self.funkcja_c(self.wektor_q, wekt_2, "+")
        tmp3 = self.lista_wektorow - wektor_2 + self.shift
        tmp4 = self.magnetyzacja[tmp3[:, 0], tmp3[:, 1]]
        return tmp1 * tmp2 * tmp4 / self.H0

    def trzecie_wyrazenie_yx(self, wektor_2):
        wekt_2 = 2 * np.pi * wektor_2 / self.a
        tmp1 = self.funkcja_c(self.wektor_q, wekt_2, "+")
        tmp3 = self.lista_wektorow - wektor_2 + self.shift
        tmp4 = self.magnetyzacja[tmp3[:, 0], tmp3[:, 1]]
        return tmp1 * tmp4 / self.H0

    # noinspection PyTypeChecker
    def czwarte_wyrazenie(self, wektor_2):
        """
        Metoda obliczająca człon elmentu macierzowego pochodzenia dipolowego.
        :param wektor_1: i-ty wektor.
        :param wektor_2: j-ty wektor.
        :return: Wynikiem jest czwarte wyrażenie w sumie na element macierzy M.
        """
        wektor_1 = self.lista_wektorow
        wekt_2 = 2 * np.pi * wektor_2 / self.a
        wekt_1 = 2 * np.pi * wektor_1 / self.a
        tmp1 = (wekt_1[:, 0] + wekt_2[0]) ** 2 / np.linalg.norm(wekt_1 - wekt_2) ** 2
        tmp2 = 1 - self.funkcja_c(wektor_1, wekt_2, "-")
        tmp3 = wektor_1 - wektor_2 + self.shift
        tmp4 = self.magnetyzacja[tmp3[:, 0], tmp3[:, 1]]
        return tmp1 * tmp2 * tmp4 / self.H0

    def wypelnienie_macierzy(self):
        """
        Główna metoda tej klasy. Wywołuje ona dwie metody: 'macierz_xy' oraz 'macierz_yx. W pętli, dla każdego elementu
        z odpowiednich macierzy blokowych wypełnia je.
        :param wektor_q: Blochowski wektor. Jest on "uciąglony". Jest on zmienną przy wyznaczaniu dyspersji.
        :return: Tablica do zagadnienia własnego.
        """
        self.delta_kroneckera()
        w1 = self.lista_wektorow

        pool = Pool()
        # tmp1 = np.array(pool.map(self.pole_wymiany_II, w1))
        tmp1 = np.apply_along_axis(self.pole_wymiany_II, 1, w1)
        tmp4 = np.array(pool.map(self.czwarte_wyrazenie, w1))
        tmp2 = np.array(pool.map(self.trzecie_wyrazenie_xy, w1))
        tmp3 = np.array(pool.map(self.trzecie_wyrazenie_yx, w1))
        self.macierz_M[self.ilosc_wektorow:, 0:self.ilosc_wektorow] += -tmp1 - tmp3 + tmp4
        self.macierz_M[0:self.ilosc_wektorow, self.ilosc_wektorow:] += tmp1 + tmp2 - tmp4
        pool.close()
        return self.macierz_M

    @do_cprofile
    def wypisz_macierz(self):
        """
        :return: Wypisuje tablice do pliku tekstowego.
         Ważne! Przed wypisaniem, należy wypełnić macierz_M metodą 'wypełnienie_macierzy'
        """
        self.wypelnienie_macierzy()
        np.savetxt('macierz.txt', np.array(self.macierz_M))


if __name__ == "__main__":
    q = MacierzDoZagadnienia('radius100.txt', 441, np.array([1e-9, 0]))
    q.wypisz_macierz()
