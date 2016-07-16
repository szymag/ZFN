# -*- coding: utf-8 -*-
import numpy as np
import multiprocessing as mp
from src.drawing.ParametryMaterialowe import ParametryMaterialowe
from src.eig_problem.FFTfromFile1D import FFTfromFile1D
from src.eig_problem.WektorySieciOdwrotnej import WektorySieciOdwrotnej


class MacierzDoZagadnienia(ParametryMaterialowe):
    """
    Klasa, w której tworzona jest macierz zagadnienia własnego. Pobiera ona współczynniki Fouriera z klasy FFTfromFile.
    Współczynniki są dla niej tworzone poprzez FFT z pliku graficznego.
    """

    def __init__(self):
        ParametryMaterialowe.__init__(self)
        self.macierz_M = np.zeros((2 * self.ilosc_wektorow, 2 * self.ilosc_wektorow), dtype=complex)
        self.tmp = FFTfromFile1D()
        self.magnetyzacja = self.tmp.fourier_coefficient(self.MoA, self.MoB)
        self.dlugosc_wymiany = self.tmp.fourier_coefficient(self.lA, self.lB)
        self.lista_wektorow = WektorySieciOdwrotnej(self.a, self.b, self.ilosc_wektorow).lista_wektorow1d('min')
        self.shift = len(self.lista_wektorow) - 1

    def funkcja_c(self, wektor_1, wektor_2):
        """
        Metoda obliczająca wartość funkcji C zdefinoweanej wzorem: f(g, x) = cosh(|g|x)*exp(-|g|d/2), gdzie g jest
        wektorem, a x współrzędną iksową, tzn. miejscem na warstwie, dla którego wykreśla się zależność dyspersyjną.
        :type wektor_1: Real
        :type wektor_2: Real
        :param wektor_1: Pierwszy wektor do obliczenia różnicy.
        :param wektor_2: Drugi wektor do obliczenia różnicy.
        :return: wartość funkcji C.
        """
        return np.exp(-abs((wektor_1 + wektor_2)) * self.d)


    def delta_kroneckera(self):
        """
        Metoda dodająca do odpowienich elementów macierzowych pierwszy element z wyrażenia na M: 1 lub -1
        """
        for i in range(self.ilosc_wektorow, 2 * self.ilosc_wektorow):
            self.macierz_M[i - self.ilosc_wektorow][i] += 1.
            self.macierz_M[i][i - self.ilosc_wektorow] -= 1.

    def pole_wymiany_II(self, wektor_1, wektor_2, wektor_q):
        """
        :param wektor_1: i-ty wektor.
        :param wektor_2: j-ty wektor.
        :param wektor_q: Blochowski wektor. Jest on "uciąglony". Jest on zmienną przy wyznaczaniu dyspersji.
        :return:
        """

        vec_l = np.array(self.lista_wektorow)
        wektor_2 = self.dlugosc_wymiany[vec_l - wektor_2 + self.shift] * \
                   (wektor_q + 6 *wektor_2 / self.a) * (6 * vec_l / self.a + wektor_q) *\
                   self.magnetyzacja[wektor_1 - vec_l + self.shift]
        return wektor_2 / self.H0

    def trzecie_wyrazenie_xy(self, wektor_1, wektor_2, wektor_q):
        """
        Metoda obliczająca człon elmentu macierzowego pochodzenia dipolowego.
        :param wektor_1: i-ty wektor.
        :param wektor_2: j-ty wektor.
        :param wektor_q: Blochowski wektor. Jest on "uciąglony". Jest on zmienną przy wyznaczaniu dyspersji.
        :return: Wynikiem jest drugi wyraz sumy.
        """

        tmp3 = self.magnetyzacja[wektor_1 - wektor_2 + self.shift]
        tmp2 = self.funkcja_c(wektor_q, (6 * wektor_2 / self.a))
        return  tmp3 * (1 - tmp2) / self.H0

    def trzecie_wyrazenie_yx(self, wektor_1, wektor_2, wektor_q):
        """
        Metoda obliczająca człon elmentu macierzowego pochodzenia dipolowego.
        :param wektor_1: i-ty wektor.
        :param wektor_2: j-ty wektor.
        :param wektor_q: Blochowski wektor. Jest on "uciąglony". Jest on zmienną przy wyznaczaniu dyspersji.
        :return: Wynikiem jest drugi wyraz sumy.
        """

        tmp3 = self.magnetyzacja[wektor_1 - wektor_2 + self.shift]
        tmp2 = self.funkcja_c(wektor_q, (6 * wektor_2 / self.a))
        return  tmp2 * tmp3 / self.H0

    def wypelnienie_macierzy(self, wektor_q):
        """
        Główna metoda tej klasy. Wywołuje ona dwie metody: 'macierz_xy' oraz 'macierz_yx. W pętli, dla każdego elementu
        z odpowiednich macierzy blokowych wypełnia je.
        :param wektor_q: Blochowski wektor. Jest on "uciąglony". Jest on zmienną przy wyznaczaniu dyspersji.
        :return: Tablica do zagadnienia własnego.
        """
        # TODO: Zaktualizować opis metody
        indeks = self.ilosc_wektorow
        assert len(self.lista_wektorow) == indeks, 'number of vector do not fit to matrix'
        self.delta_kroneckera()
        for i in range(indeks, 2 * indeks):
            w1 = self.lista_wektorow[i - indeks]
            w2 = self.lista_wektorow
            tmp1 = self.pole_wymiany_II(w1, w2, wektor_q)
            tmp2 = self.trzecie_wyrazenie_xy(w1, w2, wektor_q)
            tmp3 = self.trzecie_wyrazenie_yx(w1, w2, wektor_q)
            self.macierz_M[i][np.arange(indeks)] += -tmp1 - tmp2  # yx
            self.macierz_M[i - indeks][np.arange(indeks, 2 * indeks)] += tmp1 + tmp3  # xy
        return self.macierz_M

    def wypisz_macierz(self):
        """
        :return: Wypisuje tablice do pliku tekstowego.
         Ważne! Przed wypisaniem, należy wypełnić macierz_M metodą 'wypełnienie_macierzy'
        """
        self.wypelnienie_macierzy(1e-9)
        np.savetxt('macierz.txt', np.array(self.macierz_M))

#MacierzDoZagadnienia().wypisz_macierz()