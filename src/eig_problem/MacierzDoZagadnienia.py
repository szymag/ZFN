# -*- coding: utf-8 -*-

from math import exp, cosh, sqrt
import numpy as np
from src.eig_problem.FFTfromFile import FFTfromFile
from src.eig_problem.ParametryMaterialowe import ParametryMaterialowe
from src.eig_problem.WektorySieciOdwrotnej import WektorySieciOdwrotnej


class MacierzDoZagadnienia:
    """
    Klasa, w której tworzona jest macierz zagadnienia własnego. Pobiera ona współczynniki Fouriera z klasy FFTfromFile.
    Współczynniki są dla niej tworzone poprzez FFT z pliku graficznego.
    """

    def __init__(self, input_fft, ilosc_wektorow, a=ParametryMaterialowe.a, MoA=ParametryMaterialowe.MoA,
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

    def norma_wektorow(self, wektor_1, wektor_2, znak):
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
    def pole_wymiany_II(self, wektor_1, wektor_2, wektor_q):
        vec_l = self.lista_wektorow
        #TODO: Zwryfikować dl_wym
        dl_wym = np.sum(self.dlugosc_wymiany[vec_l - np.transpose(wektor_2) + np.transpose(self.shift)])
        skalarny = np.dot(wektor_q + 2 * np.pi * vec_l / self.a, wektor_q + 2 * np.pi * wektor_2 / self.a)
        vec_l = np.repeat(vec_l, self.ilosc_wektorow, axis=0).reshape((self.ilosc_wektorow, self.ilosc_wektorow, 2))
        wektor_1 = np.repeat(wektor_1[np.newaxis, :], self.ilosc_wektorow, 0)
        mag = np.sum(self.magnetyzacja[list(np.transpose(wektor_1 - vec_l + self.shift, axes=[2, 0, 1]))], axis=0)
        return skalarny * mag * dl_wym / self.H0

    # noinspection PyTypeChecker
    def trzecie_wyrazenie_xy(self, wektor_1, wektor_2, wektor_q):
        """
        Metoda obliczająca człon elmentu macierzowego pochodzenia dipolowego.
        :param wektor_1: i-ty wektor.
        :param wektor_2: j-ty wektor.
        :param wektor_q: Blochowski wektor. Jest on "uciąglony". Jest on zmienną przy wyzn aczaniu dyspersji.
        :return: Wynikiem jest drugi wyraz sumy.
        """
        tmp1 = np.transpose(wektor_q + 2 * np.pi * wektor_2 / self.a)[1] ** 2/ \
               np.linalg.norm(wektor_q + 2 * np.pi * wektor_2 / self.a) ** 2
        tmp2 = 1 - self.funkcja_c(wektor_q, 2 * np.pi * wektor_2 / self.a, "+")
        tmp3 = self.magnetyzacja[list(np.transpose(wektor_1 - wektor_2 + self.shift))]
        return tmp1 * tmp2 * tmp3 / self.H0

    def trzecie_wyrazenie_yx(self, wektor_1, wektor_2, wektor_q):
        tmp1 = 1 - self.funkcja_c(wektor_q, 2 * np.pi * wektor_2 / self.a, "+")
        tmp2 = self.magnetyzacja[list(np.transpose(wektor_1 - wektor_2 + self.shift))]
        return tmp1 * tmp2 / self.H0

    # noinspection PyTypeChecker
    def czwarte_wyrazenie(self, wektor_1, wektor_2):
        """
        Metoda obliczająca człon elmentu macierzowego pochodzenia dipolowego.
        :param wektor_1: i-ty wektor.
        :param wektor_2: j-ty wektor.
        :return: Wynikiem jest czwarte wyrażenie w sumie na element macierzy M.
        """
        tmp = self.magnetyzacja[list(np.transpose(wektor_1 - wektor_2 + self.shift))]
        tmp1 = 1 - self.funkcja_c(wektor_1, wektor_2, "-")
        tmp2 = np.transpose(2 * np.pi * wektor_1 / self.a + 2 * np.pi * wektor_2 / self.a)[0] ** 2 / \
               np.linalg.norm(2 * np.pi * wektor_1 / self.a - 2 * np.pi * wektor_2 / self.a) ** 2
        return tmp * tmp1 * tmp2 / self.H0

    def wypelnienie_macierzy(self, wektor_q):
        """
        Główna metoda tej klasy. Wywołuje ona dwie metody: 'macierz_xy' oraz 'macierz_yx. W pętli, dla każdego elementu
        z odpowiednich macierzy blokowych wypełnia je.
        :param wektor_q: Blochowski wektor. Jest on "uciąglony". Jest on zmienną przy wyznaczaniu dyspersji.
        :return: Tablica do zagadnienia własnego.
        """
        self.delta_kroneckera()
        w1 = self.lista_wektorow
        for i in range(self.ilosc_wektorow):
            tmp1 = self.pole_wymiany_II(w1, w1[i], wektor_q)
            tmp2 = self.trzecie_wyrazenie_xy(w1, w1[i], wektor_q)
            tmp3 = self.trzecie_wyrazenie_yx(w1, w1[i], wektor_q)
            tmp4 = self.czwarte_wyrazenie(w1, w1[i])
            self.macierz_M[i + self.ilosc_wektorow][np.arange(self.ilosc_wektorow)] += -tmp1 - tmp3 + tmp4
            self.macierz_M[i][np.arange(self.ilosc_wektorow, 2 * self.ilosc_wektorow)] += tmp1 + tmp2 - tmp4
        return self.macierz_M

    def wypisz_macierz(self):
        """
        :return: Wypisuje tablice do pliku tekstowego.
         Ważne! Przed wypisaniem, należy wypełnić macierz_M metodą 'wypełnienie_macierzy'
        """
        self.wypelnienie_macierzy([1e-9, 0])
        np.savetxt('macierz.txt', np.array(self.macierz_M))


if __name__ == "__main__":
    q = MacierzDoZagadnienia('radius100.txt', 9)
    q.wypisz_macierz()