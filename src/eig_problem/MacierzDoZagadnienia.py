# -*- coding: utf-8 -*-
from ctypes import *
from math import exp, cosh, sqrt
import numpy as np
from src.eig_problem.cProfiler import do_cprofile
from src.eig_problem.DFT import DFT
from src.eig_problem.FFTfromFile import FFTfromFile
from src.eig_problem.ParametryMaterialowe import ParametryMaterialowe
from src.eig_problem.WektorySieciOdwrotnej import WektorySieciOdwrotnej
import multiprocessing as mp
class vec2d(Structure):
    _fields_ = [("x_", c_longlong),
                ("y_", c_longlong)]

class vec2d_float(Structure):
    _fields_ = [("x_", c_double),
                ("y_", c_double)]



class MacierzDoZagadnienia(ParametryMaterialowe):
    """
    Klasa, w której tworzona jest macierz zagadnienia własnego. Pobiera ona współczynniki Fouriera z klasy FFTfromFile.
    Współczynniki są dla niej tworzone poprzez FFT z pliku graficznego.
    """

    def __init__(self, ilosc_wektorow, skad_wspolczynnik, typ_pola_wymiany):
        ParametryMaterialowe.__init__(self, ilosc_wektorow, typ_pola_wymiany)
        self.macierz_M = np.zeros((2 * ilosc_wektorow, 2 * ilosc_wektorow), dtype=complex)
        self.ilosc_wektorow = ilosc_wektorow
        self.typ_pola_wymiany = typ_pola_wymiany
        if skad_wspolczynnik == 'FFT':
            self.tmp = FFTfromFile(ilosc_wektorow, typ_pola_wymiany)
            self.slownik_magnetyzacja = self.tmp.fourier_coefficient()
            self.slownik_dlugosc_wymiany = self.tmp.exchange_length()
        else:
            wspolczynniki = DFT(self.ilosc_wektorow, typ_pola_wymiany).slownik_wspolczynnikow()
            self.slownik_dlugosc_wymiany = wspolczynniki[1]
            self.slownik_magnetyzacja = wspolczynniki[0]
        self.lista_wektorow = WektorySieciOdwrotnej(self.a, self.b, ilosc_wektorow).lista_wektorow('min')

        self.pole_wymiany_dll = CDLL('./pole_wymiany_c/pole_wymiany_ii.so')
        for k in self.slownik_dlugosc_wymiany.keys():
            v = self.slownik_dlugosc_wymiany[k]
            self.pole_wymiany_dll.add_dlugosc_wymiany_value(c_longlong(k[0]),
                    c_longlong(k[1]), c_double(v.real), c_double(v.imag))
        for k in self.slownik_magnetyzacja.keys():
            v = self.slownik_magnetyzacja[k]
            self.pole_wymiany_dll.add_magnetyzacja_value(c_longlong(k[0]),
                    c_longlong(k[1]), c_double(v.real), c_double(v.imag))

        vec2d_array_type = vec2d * len(self.lista_wektorow)
        arr = vec2d_array_type(*self.lista_wektorow)
        self.pole_wymiany_dll.init_lista_wektorow(arr, len(self.lista_wektorow), c_double(self.H0))
        self.pole_wymiany_dll.tmp_value.argtypes = [POINTER(vec2d),
                                                    POINTER(vec2d),
                                                    POINTER(vec2d_float),
                                                    POINTER(vec2d_float)]
        self.pole_wymiany_dll.tmp_value.restype = None

    def magnetyzacja(self, wektor_1, wektor_2):
        """
        Metoda wyliczająca współczynnik Fouriera. Jako argumenty podawne są dwa wektory i wyliczana jest
        z nich różnica.
        :type wektor_1: tuple
        :type wektor_2: tuple
        :param wektor_1: Pierwszy wektor do obliczenia różnicy.
        :param wektor_2: Drugi wektor do obliczenia różnicy.
        :return: współczynnik Fouriera dla różnicy wektorów sieci odwrotnej.
        """
        return self.slownik_magnetyzacja[(wektor_1[0] - wektor_2[0], wektor_1[1] - wektor_2[1])]

    def dlugosc_wymiany(self, wektor_1, wektor_2):
        """
        Metoda obliczająca długość wymiany, dla dwóch zadanych wektorów.
        z nich różnica.
        :type wektor_1: tuple
        :type wektor_2: tuple
        :param wektor_1: Pierwszy wektor do obliczenia różnicy.
        :param wektor_2: Drugi wektor do obliczenia różnicy.
        :return: Długość wymiany w postaci odpowiadającego różnicy wektorów współczynnika Fouriera.
        """
        return self.slownik_dlugosc_wymiany[wektor_1[0] - wektor_2[0], wektor_1[1] - wektor_2[1]]

    def funkcja_c(self, wektor_1, wektor_2, znak):
        """
        Metoda obliczająca wartość funkcji C zdefinoweanej wzorem: f(g, x) = cosh(|g|x)*exp(-|g|d/2), gdzie g jest
        wektorem, a x współrzędną iksową, tzn. miejscem na warstwie, dla którego wykreśla się zależność dyspersyjną.
        :type wektor_1: tuple
        :type wektor_2: tuple
        :type znak: str
        :param wektor_1: Pierwszy wektor do obliczenia różnicy.
        :param wektor_2: Drugi wektor do obliczenia różnicy.
        :param znak: określa czy obliczana ma być różnica, czy suma wektorów.
        :return: wartość funkcji C.
        """
        wekt_wypadkowy = self.norma_wektorow(wektor_1, wektor_2, znak)
        return cosh(wekt_wypadkowy * self.x) * exp(-wekt_wypadkowy * self.d / 2.)

    def norma_wektorow(self, wektor_1, wektor_2, znak):
        """
        Metoda wyliczająca wypadkową długość dwóch wektorów.
        :type wektor_1: tuple
        :type wektor_2: tuple
        :type znak: str
        :param wektor_1: Pierwszy wektor do obliczenia różnicy.
        :param wektor_2: Drugi wektor do obliczenia różnicy.
        :param znak: Określa czy obliczana ma być różnica, czy suma wektorów.
        :return: W zlależności od znaku, zwraca normę z sumy, lub różnicy wektorów.
        """
        if znak == "-":
            return sqrt((wektor_1[0] - wektor_2[0])**2 + (wektor_1[1] - wektor_2[1])**2)
        elif znak == "+":
            return sqrt((wektor_1[0] + wektor_2[0])**2 + (wektor_1[1] + wektor_2[1])**2)

    def delta_kroneckera(self):
        """
        Metoda dodająca do odpowienich elementów macierzowych pierwszy element z wyrażenia na M: 1 lub -1
        """
        for i in range(self.ilosc_wektorow, 2 * self.ilosc_wektorow):
            self.macierz_M[i - self.ilosc_wektorow][i] += 1.
            self.macierz_M[i][i - self.ilosc_wektorow] -= 1.

    def pole_wymiany_I(self, wektor_1, wektor_2, wektor_q):
        """
        Metoda obliczająca człon elmentu macierzowego związana z polem wymiany.
        :type wektor_1: tuple
        :type wektor_2: tuple
        :type wektor_q: tuple
        :param wektor_1: i-ty wektor
        :param wektor_2: j-ty wektor
        :param wektor_q: Blochowski wektor. Jest on "uciąglony". Jest on zmienną przy wyznaczaniu dyspersji.
        :return: Wynikiem jest drugi wyraz sumy.
        """
        tmp1 = (wektor_q[0] + wektor_2[0]) * (wektor_q[0] + wektor_1[0]) + \
               (wektor_q[1] + wektor_1[1]) * (wektor_q[1] + wektor_2[1])
        tmp2 = self.dlugosc_wymiany(wektor_1, wektor_2)

        return tmp1 * tmp2 / self.H0

    def pole_wymiany_II(self, wektor_1, wektor_2, wektor_q):
        res_tmp = vec2d_float()
        self.pole_wymiany_dll.tmp_value(byref(vec2d(*wektor_1)),
                                        byref(vec2d(*wektor_2)),
                                        byref(vec2d_float(*wektor_q)),
                                        byref(res_tmp))
        return complex(res_tmp.x_, res_tmp.y_)

    def trzecie_wyrazenie_xy_yx(self, wektor_1, wektor_2, wektor_q):
        """
        Metoda obliczająca człon elmentu macierzowego pochodzenia dipolowego.
        :type wektor_1: tuple
        :type wektor_2: tuple
        :type wektor_q: tuple
        :type typ_macierzy: str
        :param wektor_1: i-ty wektor.
        :param wektor_2: j-ty wektor.
        :param wektor_q: Blochowski wektor. Jest on "uciąglony". Jest on zmienną przy wyznaczaniu dyspersji.
        :param typ_macierzy: Określa, do której z macierzy blokowych odnosi się wyrażenie.
        :return: Wynikiem jest drugi wyraz sumy.
        """
        tmp1 = self.norma_wektorow(wektor_q, wektor_2, '+')
        assert tmp1 != 0., 'probably insert forbidden q vector e.g. q = 0, q = 1'
        tmp2 = self.funkcja_c(wektor_q, wektor_2, "+")
        tmp3 = self.magnetyzacja(wektor_1, wektor_2)
        return (wektor_q[0] + wektor_2[0]) ** 2 / (self.H0 * tmp1 ** 2) * (1 - tmp2) * tmp3, tmp2 * tmp3 / self.H0

    def czwarte_wyrazenie(self, wektor_1, wektor_2):
        """
        Metoda obliczająca człon elmentu macierzowego pochodzenia dipolowego.
        :type wektor_1: tuple
        :type wektor_2: tuple
        :param wektor_1: i-ty wektor.
        :param wektor_2: j-ty wektor.
        :return: Wynikiem jest czwarte wyrażenie w sumie na element macierzy M.
        """
        return (wektor_1[1] - wektor_2[1]) ** 2 / \
               (1e-36 + self.H0 * self.norma_wektorow(wektor_1, wektor_2, "-") ** 2) * \
               self.magnetyzacja(wektor_1, wektor_2) * (1 - self.funkcja_c(wektor_1, wektor_2, "-"))


    def wypelnienie_macierzy(self, wektor_q):
        """
        Główna metoda tej klasy. Wywołuje ona dwie metody: 'macierz_xy' oraz 'macierz_yx. W pętli, dla każdego elementu
        z odpowiednich macierzy blokowych wypełnia je.
        :param wektor_q: Blochowski wektor. Jest on "uciąglony". Jest on zmienną przy wyznaczaniu dyspersji.
        :return: Tablica do zagadnienia własnego.
        """
        # TODO: Zaktualizować opis metody
        assert type(wektor_q) == tuple, \
            'form of wektor_q is forbidden. wektor_q should be touple'
        assert len(wektor_q) == 2, \
            'form of wektor_q is forbidden. wektor_q should have two arguments'
        indeks = self.ilosc_wektorow
        assert len(self.lista_wektorow) == indeks, 'number of vector do not fit to matrix'
        self.delta_kroneckera()

        for i in range(indeks, 2 * indeks):

            w1 = self.lista_wektorow[i - indeks]

            for j in range(0, indeks):
                w2 = self.lista_wektorow[j]
                tmp1 = self.pole_wymiany_II(w1, w2, wektor_q)
                tmp2 = self.trzecie_wyrazenie_xy_yx(w1, w2, wektor_q)
                tmp4 = self.czwarte_wyrazenie(w1, w2)

                self.macierz_M[i][j] += -tmp1 - tmp2[1] + tmp4
                self.macierz_M[i - indeks][j + indeks] += tmp1 + tmp2[0] - tmp4
        return self.macierz_M

    def wypisz_macierz(self):
        """
        :return: Wypisuje tablice do pliku tekstowego.
         Ważne! Przed wypisaniem, należy wypełnić macierz_M metodą 'wypełnienie_macierzy'
        """
        np.savetxt('macierz.txt', np.array(self.macierz_M))
