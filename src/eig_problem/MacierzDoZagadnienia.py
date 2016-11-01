# -*- coding: utf-8 -*-
import numpy as np

from src.eig_problem.FFTfromFile1D import FFTfromFile1D
from src.eig_problem.WektorySieciOdwrotnej import WektorySieciOdwrotnej
from src.eig_problem.ParametryMaterialowe import ParametryMaterialowe

class MacierzDoZagadnienia:
    """
    Klasa, w której tworzona jest macierz zagadnienia własnego. Pobiera ona współczynniki Fouriera z klasy FFTfromFile.
    Współczynniki są dla niej tworzone poprzez FFT z pliku graficznego.
    """

    def __init__(self, input_fft, wektor_q, a=ParametryMaterialowe.a, MoA=ParametryMaterialowe.MoA,
                 MoB=ParametryMaterialowe.MoB, lA=ParametryMaterialowe.lA,
                 lB=ParametryMaterialowe.lB, d=ParametryMaterialowe.d,
                x=ParametryMaterialowe.x, H0=ParametryMaterialowe.H0):

        self.a = a
        self.d = d
        self.x = x
        self.H0 = H0
        self.tmp = FFTfromFile1D(input_fft)
        self.ilosc_wektorow = self.tmp.ilosc_wektorow
        self.macierz_M = np.zeros((2 * self.ilosc_wektorow, 2 * self.ilosc_wektorow), dtype=complex)
        self.magnetyzacja = self.tmp.fourier_coefficient(MoA, MoB)
        self.dlugosc_wymiany = self.tmp.fourier_coefficient(lA, lB)
        self.lista_wektorow = WektorySieciOdwrotnej(self.ilosc_wektorow).lista_wektorow1d('min')
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
        return np.exp(-abs((wektor_1 + wektor_2)) * self.d/2)

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
        vec_l = np.transpose(np.broadcast_to(self.lista_wektorow , (self.ilosc_wektorow, self.ilosc_wektorow)))
        tmp1 = self.dlugosc_wymiany[vec_l - wektor_2 + self.shift]
        tmp2 = self.magnetyzacja[wektor_1 - vec_l + self.shift]
        tmp3 = (wektor_q + 2 * np.pi * wektor_2 / self.a) * (2 * np.pi * vec_l / self.a + wektor_q)
        return np.sum(tmp1 * tmp2 * tmp3, axis=0) / self.H0

    def trzecie_wyrazenie_xy(self, wektor_1, wektor_2, wektor_q):
        """
        Metoda obliczająca człon elmentu macierzowego pochodzenia dipolowego.
        :param wektor_1: i-ty wektor.
        :param wektor_2: j-ty wektor.
        :param wektor_q: Blochowski wektor. Jest on "uciąglony". Jest on zmienną przy wyznaczaniu dyspersji.
        :return: Wynikiem jest drugi wyraz sumy.
        """

        tmp3 = self.magnetyzacja[wektor_1 - wektor_2 + self.shift]
        tmp2 = 1 - self.funkcja_c(wektor_q, (2 * np.pi * wektor_2 / self.a))
        return tmp3 * tmp2 / self.H0

    def trzecie_wyrazenie_yx(self, wektor_1, wektor_2, wektor_q):
        """
        Metoda obliczająca człon elmentu macierzowego pochodzenia dipolowego.
        :param wektor_1: i-ty wektor.
        :param wektor_2: j-ty wektor.
        :param wektor_q: Blochowski wektor. Jest on "uciąglony". Jest on zmienną przy wyznaczaniu dyspersji.
        :return: Wynikiem jest drugi wyraz sumy.
        """
        #print(wektor_1)
        tmp3 = self.magnetyzacja[wektor_1 - wektor_2 + self.shift]
        tmp2 = self.funkcja_c(wektor_q, (2 * np.pi * wektor_2 / self.a))
        return tmp2 * tmp3 / self.H0

    def wypelnienie_macierzy(self, wektor_q):
        """
        Główna metoda tej klasy. Wywołuje ona dwie metody: 'macierz_xy' oraz 'macierz_yx. W pętli, dla każdego elementu
        z odpowiednich macierzy blokowych wypełnia je.
        :param wektor_q: Blochowski wektor. Jest on "uciąglony". Jest on zmienną przy wyznaczaniu dyspersji.
        :return: Tablica do zagadnienia własnego.
        """
        # TODO: Zaktualizować opis metody
        indeks = self.ilosc_wektorow
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

if __name__ == "__main__":
    q = MacierzDoZagadnienia('p_coef_10*2.txt', 1e-9)
    q.wypisz_macierz()