# -*- coding: utf-8 -*-

from math import exp, cosh, sqrt

from numpy import dot
from numpy import zeros, array, savetxt

from src.eig_problem.DFT import DFT
from src.eig_problem.FFTfromFile import FFTfromFile
from src.eig_problem.ParametryMaterialowe import ParametryMaterialowe
from src.eig_problem.WektorySieciOdwrotnej import WektorySieciOdwrotnej


class MacierzDoZagadnienia(ParametryMaterialowe):
    """
    Klasa, w której tworzona jest macierz zagadnienia własnego. Pobiera ona współczynniki Fouriera z klasy FFTfromFile.
    Współczynniki są dla niej tworzone poprzez FFT z pliku graficznego.
    """

    def __init__(self, ilosc_wektorow, skad_wspolczynnik, typ_pola_wymiany):
        ParametryMaterialowe.__init__(self, ilosc_wektorow, typ_pola_wymiany)
        self.macierz_M = zeros((2 * ilosc_wektorow, 2 * ilosc_wektorow), dtype=complex)
        self.ilosc_wektorow = ilosc_wektorow
        self.typ_pola_wymiany = typ_pola_wymiany
        if skad_wspolczynnik == 'FFT':
            self.slownik_wspolczynnik = FFTfromFile(ilosc_wektorow, typ_pola_wymiany).fourier_coefficient()
            self.slownik_dlugosc_wymiany = FFTfromFile(ilosc_wektorow, typ_pola_wymiany).exchange_length()
            self.lista_wektorow = FFTfromFile(ilosc_wektorow, typ_pola_wymiany).vector_to_matrix()
        else:
            self.slownik_dlugosc_wymiany = DFT(self.ilosc_wektorow, typ_pola_wymiany).slownik_wspolczynnikow()[1]
            self.slownik_wspolczynnik = DFT(self.ilosc_wektorow, typ_pola_wymiany).slownik_wspolczynnikow()[0]
            self.lista_wektorow = WektorySieciOdwrotnej(self.a, self.b, ilosc_wektorow).lista_wektorow('min')

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
        return self.slownik_wspolczynnik[(wektor_1[0] - wektor_2[0], wektor_1[1] - wektor_2[1])]

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

    @staticmethod
    def suma_roznica_wektorow(wektor_1, wektor_2, znak):
        """
        Metoda, która w zależności od znaku oblicza sumę, bądż różnicę wektorów.
        :type wektor_1: tuple
        :type wektor_2: tuple
        :type znak: str
        :param wektor_1: Pierwszy wektor do obliczenia różnicy.
        :param wektor_2: Drugi wektor do obliczenia różnicy.
        :param znak: Określa czy obliczana ma być różnica, czy suma wektorów.
        :return: suma lub rożnica wektorów
        """
        assert type(znak) == str, \
            'znak is sign between two vector. Should be string'
        assert znak == '+' or znak == '-', \
            'only - and + are permitted'

        if znak == "-":
            return (wektor_1[0] - wektor_2[0], wektor_1[1] - wektor_2[1])
        elif znak == "+":
            return (wektor_1[0] + wektor_2[0], wektor_1[1] + wektor_2[1])

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
        assert type(znak) == str, \
            'znak is sign between two vector. Should be string'
        assert znak == '+' or znak == '-', \
            'only - and + are permitted'

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
        assert type(znak) == str, \
            'znak is sign between two vector. Should be string'
        assert znak == '+' or znak == '-', \
            'only - and + are permitted'
        wekt_wypadkowy = self.suma_roznica_wektorow(wektor_1, wektor_2, znak)
        return sqrt(wekt_wypadkowy[0] ** 2 + wekt_wypadkowy[1] ** 2)

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
        tmp1 = (wektor_q[0] + wektor_2[0]) * (wektor_q[0] + wektor_1[0]) \
               + (wektor_q[1] + wektor_1[1]) * (wektor_q[1] + wektor_2[1])
        tmp2 = self.dlugosc_wymiany(wektor_1, wektor_2)
        return tmp1 * tmp2 / self.H0

    def pole_wymiany_II(self, wektor_1, wektor_2, wektor_q):
        """
        :type wektor_1: tuple
        :type wektor_2: tuple
        :type wektor_q: tuple
        :param wektor_1: i-ty wektor.
        :param wektor_2: j-ty wektor.
        :param wektor_q: Blochowski wektor. Jest on "uciąglony". Jest on zmienną przy wyznaczaniu dyspersji.
        :return:
        """
        # TODO: Usprawnić iloczyn skalarny
        tmp = 0.
        for vec_l in self.lista_wektorow:
            tmp += dot((wektor_q[0] + wektor_2[0], wektor_q[1] + wektor_2[1]),
                       (wektor_q[0] + vec_l[0], wektor_q[1] + vec_l[1])) * \
                   self.dlugosc_wymiany(vec_l, wektor_2) * self.magnetyzacja(wektor_1, vec_l) / self.H0
        return tmp

    def trzecie_wyrazenie(self, wektor_1, wektor_2, wektor_q, typ_macierzy):
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
        assert typ_macierzy == 'xy' or typ_macierzy == 'yx', \
            'it is assumed that block matrixes are named xy or yx'
        tmp1 = self.norma_wektorow(wektor_q, wektor_2, '+')
        assert tmp1 != 0., 'probably insert forbidden q vector e.g. q = 0, q = 1'
        tmp2 = self.funkcja_c(wektor_q, wektor_2, "+")
        tmp3 = self.magnetyzacja(wektor_1, wektor_2)
        if typ_macierzy == 'xy':
            return (wektor_q[0] + wektor_2[0]) ** 2 / (self.H0 * tmp1 ** 2) * (1 - tmp2) * tmp3
        elif typ_macierzy == 'yx':
            return tmp2 * tmp3 / self.H0

    def czwarte_wyrazenie(self, wektor_1, wektor_2):
        """
        Metoda obliczająca człon elmentu macierzowego pochodzenia dipolowego.
        :type wektor_1: tuple
        :type wektor_2: tuple
        :param wektor_1: i-ty wektor.
        :param wektor_2: j-ty wektor.
        :return: Wynikiem jest czwarte wyrażenie w sumie na element macierzy M.
        """
        return (wektor_1[1] - wektor_2[1]) ** 2 / (1e-36 + self.H0 * self.norma_wektorow(wektor_1, wektor_2, "-") ** 2) \
               * self.magnetyzacja(wektor_1, wektor_2) * (1 - self.funkcja_c(wektor_1, wektor_2, "-"))

    def macierz_xy(self, wektor_1, wektor_2, wektor_q):
        """
        Wywołuje metody wyliczające wyrażenia wchodzące do elementów macierzowych. Ta metoda zapełnia górną prawą
        macierz blokową, wchodzącą do zagadnienia własnego.
        :type wektor_1: tuple
        :type wektor_2: tuple
        :type wektor_q: tuple
        :param wektor_1: i-ty wektor.
        :param wektor_2: j-ty wektor.
        :param wektor_q: Blochowski wektor. Jest on "uciąglony". Jest on zmienną przy wyznaczaniu dyspersji.
        :return: element macierzowy
        """
        assert len(wektor_1) == 2, \
            'form of wektor_q is forbidden. wektor_1 should have two arguments'
        assert len(wektor_2) == 2, \
            'form of wektor_q is forbidden. wektor_2 should have two arguments'
        assert len(wektor_q) == 2, \
            'form of wektor_q is forbidden. wektor_q should have two arguments'
        return \
            self.pole_wymiany_I(wektor_1, wektor_2, wektor_q) \
            + self.trzecie_wyrazenie(wektor_1, wektor_2, wektor_q, "xy") \
            - self.czwarte_wyrazenie(wektor_1, wektor_2)

    def macierz_yx(self, wektor_1, wektor_2, wektor_q):
        """
        Wywołuje metody wyliczające wyrażenia wchodzące do elementów macierzowych. Ta metoda zapełnia lewą dolną
        macierz blokową, wchodzącą do zagadnienia własnego.
        :type wektor_1: tuple
        :type wektor_2: tuple
        :type wektor_q: tuple
        :param wektor_1: i-ty wektor.
        :param wektor_2: j-ty wektor.
        :param wektor_q: Blochowski wektor. Jest on "uciąglony". Jest on zmienną przy wyznaczaniu dyspersji.
        :return: element macierzowy
        """
        assert len(wektor_1) == 2, \
            'form of wektor_q is forbidden. wektor_1 should have two arguments'
        assert len(wektor_2) == 2, \
            'form of wektor_q is forbidden. wektor_2 should have two arguments'
        assert len(wektor_q) == 2, \
            'form of wektor_q is forbidden. wektor_q should have two arguments'
        return \
            - self.pole_wymiany_I(wektor_1, wektor_2, wektor_q) \
            - self.trzecie_wyrazenie(wektor_1, wektor_2, wektor_q, "yx") \
            + self.czwarte_wyrazenie(wektor_1, wektor_2)

    def wypelnienie_macierzy(self, wektor_q):
        """
        Główna metoda tej klasy. Wywołuje ona dwie metody: 'macierz_xy' oraz 'macierz_yx. W pętli, dla każdego elementu
        z odpowiednich macierzy blokowych wypełnia je.
        :param wektor_q: Blochowski wektor. Jest on "uciąglony". Jest on zmienną przy wyznaczaniu dyspersji.
        :return: Tablica do zagadnienia własnego.
        """
        assert type(wektor_q) == tuple, \
            'form of wektor_q is forbidden. wektor_q should be touple'
        assert len(wektor_q) == 2, \
            'form of wektor_q is forbidden. wektor_q should have two arguments'

        indeks = self.ilosc_wektorow
        lista_wektorow = self.lista_wektorow
        assert len(lista_wektorow) == indeks, 'number of vector do not fit to matrix'
        self.delta_kroneckera()
        for i in range(indeks, 2 * indeks):
            for j in range(0, indeks):
                w1 = lista_wektorow[i - indeks]
                w2 = lista_wektorow[j]
                tmp1 = self.pole_wymiany_I(w1, w2, wektor_q)
                tmp2 = self.trzecie_wyrazenie(w1, w2, wektor_q, "yx")
                tmp3 = self.trzecie_wyrazenie(w1, w2, wektor_q, "xy")
                tmp4 = self.czwarte_wyrazenie(w1, w2)

                self.macierz_M[i][j] += -tmp1 - tmp3 + tmp4
                self.macierz_M[i - indeks][j + indeks] += tmp1 + tmp2 - tmp4
        return self.macierz_M

    def wypisz_macierz(self):
        """
        :return: Wypisuje tablice do pliku tekstowego.
         Ważne! Przed wypisaniem, należy wypełnić macierz_M metodą 'wypełnienie_macierzy'
        """
        savetxt('macierz.txt', array(self.macierz_M))
