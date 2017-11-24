# -*- coding: utf-8 -*-
import numpy as np
from math import radians
from src.eig_problem.FFTfromFile1D import FFTfromFile1D
from src.eig_problem.WektorySieciOdwrotnej import WektorySieciOdwrotnej
from src.eig_problem.ParametryMaterialowe import ParametryMaterialowe

class MacierzDoZagadnienia:
    """
    Klasa, w której tworzona jest macierz zagadnienia własnego. Pobiera ona współczynniki Fouriera z klasy FFTfromFile.
    Współczynniki są dla niej tworzone poprzez FFT z pliku graficznego.
    """

    def __init__(self, input_fft, wektor_q, H0=ParametryMaterialowe.H0, a=ParametryMaterialowe.a,
                 MoA=ParametryMaterialowe.MoA,
                 MoB=ParametryMaterialowe.MoB, lA=ParametryMaterialowe.lA,
                 lB=ParametryMaterialowe.lB, d=ParametryMaterialowe.d,
                x=ParametryMaterialowe.x, angle=ParametryMaterialowe.angle):

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
        self.angle = angle

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
        return np.cosh(abs((wektor_1 + wektor_2)) * self.x)*np.exp(-abs((wektor_1 + wektor_2)) * self.d / 2)

    def delta_kroneckera(self):
        """
        Metoda dodająca do odpowienich elementów macierzowych pierwszy element z wyrażenia na M: 1 lub -1
        """
        for i in range(self.ilosc_wektorow, 2 * self.ilosc_wektorow):
            self.macierz_M[i - self.ilosc_wektorow][i] += 1.
            self.macierz_M[i][i - self.ilosc_wektorow] -= 1.

    def exchange_field(self, wektor_1, wektor_2, wektor_q):
        """
        :param wektor_1: i-ty wektor.
        :param wektor_2: j-ty wektor.
        :param wektor_q: Blochowski wektor. Jest on "uciąglony". Jest on zmienną przy wyznaczaniu dyspersji.
        :return:
        """
        vec_l = np.transpose(np.broadcast_to(self.lista_wektorow ,
         (self.ilosc_wektorow, self.ilosc_wektorow)))
        tmp1 = self.dlugosc_wymiany[vec_l - wektor_2 + self.shift]
        tmp2 = self.magnetyzacja[wektor_1 - vec_l + self.shift]
        tmp3 = (wektor_q + 2 * np.pi * wektor_2 / self.a) * (2 * np.pi * vec_l / self.a + wektor_q)
        return np.sum(tmp1 * tmp2 * tmp3, axis=0) / self.H0

    def dynamic_demagnetizing_field_in_plane(self, wektor_1, wektor_2, wektor_q):
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

    def dynamic_demagnetizing_field_out_of_plane(self, wektor_1, wektor_2, wektor_q):
        """
        Metoda obliczająca człon elmentu macierzowego pochodzenia dipolowego.
        :param wektor_1: i-ty wektor.
        :param wektor_2: j-ty wektor.
        :param wektor_q: Blochowski wektor. Jest on "uciąglony". Jest on zmienną przy wyznaczaniu dyspersji.
        :return: Wynikiem jest drugi wyraz sumy.
        """
        tmp3 = self.magnetyzacja[wektor_1 - wektor_2 + self.shift]
        #tmp3[self.shift] = 0
        tmp2 = self.funkcja_c(wektor_q, (2 * np.pi * wektor_2 / self.a))
        return tmp2 * tmp3 / self.H0

    def static_demagnetizing_field(self, wektor_1, wektor_2):
        """
        Metoda obliczająca człon elmentu macierzowego pochodzenia dipolowego.
        :param wektor_2: j-ty wektor.
        :return: Wynikiem jest czwarte wyrażenie w sumie na element macierzy M.
        """
        tmp1 = self.magnetyzacja[wektor_1 - wektor_2 + self.shift]
        co = np.cosh(abs((2 * np.pi * wektor_1 / self.a - 2 * np.pi * wektor_2 / self.a)) * self.x)
        tmp2 = 1 - co*np.exp(-abs((2 * np.pi * wektor_1 / self.a - 2 * np.pi * wektor_2 / self.a)) * self.d / 2)
        tmp1[self.shift] = 0
        return tmp2 * tmp1 / self.H0

    def matrix_angle_dependence(self, wektor_q):
        indeks = self.ilosc_wektorow
        self.delta_kroneckera()
        for i in range(indeks, 2 * indeks):
            w1 = self.lista_wektorow[i - indeks]
            w2 = self.lista_wektorow
            ex = self.exchange_field(w1, w2, wektor_q)
            dyn_in_plane = self.dynamic_demagnetizing_field_in_plane(w1, w2, wektor_q)
            dyn_out_plane = self.dynamic_demagnetizing_field_out_of_plane(w1, w2, wektor_q)
            static = self.static_demagnetizing_field(w1, w2)
            self.macierz_M[i][np.arange(indeks)] += -ex - dyn_in_plane * np.sin\
                (radians(self.angle))**2 + static * np.cos(radians(self.angle))**2  # yx
            self.macierz_M[i - indeks][np.arange(indeks, 2 * indeks)] += ex + dyn_out_plane\
                                                                         - static * np.cos(radians(self.angle))**2  # xy
        return self.macierz_M

    def wypisz_macierz(self):
        """
        :return: Wypisuje tablice do pliku tekstowego.
         Ważne! Przed wypisaniem, należy wypełnić macierz_M metodą 'wypełnienie_macierzy'
        """
        self.matrix_angle_dependence(1e-9)
        np.savetxt('macierz.txt', np.array(self.macierz_M))


if __name__ == "__main__":
    q = MacierzDoZagadnienia('p_coef_100*2.txt', 1e-9)
    q.wypisz_macierz()
