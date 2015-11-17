# -*- coding: utf-8 -*-


import scipy.special
from numpy import pi, sqrt, cosh, exp, dot
from numpy import zeros, array, savetxt

from src.eig_problem.ParametryMaterialowe import ParametryMaterialowe
from src.eig_problem.WektorySieciOdwrotnej import WektorySieciOdwrotnej


class MacierzDoZagadnienia(ParametryMaterialowe):
    def __init__(self, rozmiar_macierzy_blok):
        ParametryMaterialowe.__init__(self, rozmiar_macierzy_blok)
        self.macierz_M = zeros((2 * rozmiar_macierzy_blok, 2 * rozmiar_macierzy_blok))

    def wspolczynnik(self, wektor_1, wektor_2):
        """
        Metoda wyliczająca współczynnik Fouriera. Jako argumenty podawne są dwa wektory i wyliczana jest
        z nich różnica.
        :param wektor_1: Pierwszy wektor do obliczenia różnicy.
        :param wektor_2: Drugi wektor do obliczenia różnicy.
        :return: współczynnik Fouriera.
        """
        assert type(wektor_1) == tuple, \
            'form of wektor_q is forbidden. wektor_1 should be touple'
        assert len(wektor_1) == 2,\
            'form of wektor_q is forbidden. wektor_1 should have two arguments'
        assert type(wektor_2) == tuple, \
            'form of wektor_q is forbidden. wektor_2 should be touple'
        assert len(wektor_2) == 2,\
            'form of wektor_q is forbidden. wektor_2 should have two arguments'
        wekt_wypadkowy = self.suma_roznica_wektorow(wektor_1, wektor_2, '-')

        if wekt_wypadkowy[0] == 0 and wekt_wypadkowy[1] == 0:
            return (self.MoCo - self.MoPy) * pi * self.r ** 2 / (self.a ** 2) + self.MoPy
        else:
            assert wekt_wypadkowy[0] ** 2 + wekt_wypadkowy[1] ** 2 != 0, 'division by 0'
            return 2 * (self.MoCo - self.MoPy) * pi * self.r ** 2 / self.a ** 2 * \
                   scipy.special.j1(sqrt(wekt_wypadkowy[0] ** 2 + wekt_wypadkowy[1] ** 2) * self.r) \
                   / (sqrt(wekt_wypadkowy[0] ** 2 + wekt_wypadkowy[1] ** 2) * self.r)

    def wektor_pozycji(self, wektor_1, wektor_2):
        """
        Metoda wyliczająca współczynnik Fouriera. Jako argumenty podawne są dwa wektory i wyliczana jest
        z nich różnica.
        :param wektor_1: Pierwszy wektor do obliczenia różnicy.
        :param wektor_2: Drugi wektor do obliczenia różnicy.
        :return: współczynnik Fouriera.
        """
        assert type(wektor_1) == tuple, \
            'form of wektor_q is forbidden. wektor_1 should be touple'
        assert len(wektor_1) == 2,\
            'form of wektor_q is forbidden. wektor_1 should have two arguments'
        assert type(wektor_2) == tuple, \
            'form of wektor_q is forbidden. wektor_2 should be touple'
        assert len(wektor_2) == 2,\
            'form of wektor_q is forbidden. wektor_2 should have two arguments'

        zipped = list(zip(wektor_1, wektor_2))
        wekt_wypadkowy = [k[0] + k[1] for k in zipped]

        if wekt_wypadkowy[0] == 0 and wekt_wypadkowy[1] == 0:
            return (self.ACo - self.APy) * pi * self.r ** 2 / self.a ** 2 + self.APy
        else:
            return 2 * (self.ACo - self.APy) * pi * self.r ** 2 / self.a ** 2 * \
                   scipy.special.j1(sqrt(wekt_wypadkowy[0] ** 2 + wekt_wypadkowy[1] ** 2) * self.r) \
                   / (sqrt(wekt_wypadkowy[0] ** 2 + wekt_wypadkowy[1] ** 2 + (10 ** -10)) * self.r)

    @staticmethod
    def suma_roznica_wektorow(wektor_1, wektor_2, znak):
        """
        :param wektor_1: Pierwszy wektor do obliczenia różnicy.
        :param wektor_2: Drugi wektor do obliczenia różnicy.
        :param znak: Określa czy obliczana ma być różnica, czy suma wektorów.
        :return: suma lub rożnica wektorów
        """
        assert type(wektor_1) == tuple, \
            'form of wektor_q is forbidden. wektor_1 should be touple'
        assert len(wektor_1) == 2,\
            'form of wektor_q is forbidden. wektor_1 should have two arguments'
        assert type(wektor_2) == tuple, \
            'form of wektor_q is forbidden. wektor_2 should be touple'
        assert len(wektor_2) == 2,\
            'form of wektor_q is forbidden. wektor_2 should have two arguments'
        assert type(znak) == str, \
            'znak is sign between two vector. Should be string'
        assert znak == '+' or znak == '-', \
            'only - and + are permitted'

        zipped = list(zip(wektor_1, wektor_2))
        if znak == "-":
            return tuple([k[0] - k[1] for k in zipped])
        elif znak == "+":
            return tuple([k[0] + k[1] for k in zipped])

    def funkcja_c(self, wektor_1, wektor_2, znak):
        """
        Metoda obliczająca wartość funkcji C dla rożnicy wektorów
        :param wektor_1: Pierwszy wektor do obliczenia różnicy.
        :param wektor_2: Drugi wektor do obliczenia różnicy.
        :param znak: określa czy obliczana ma być różnica, czy suma wektorów.
        :return: wartość funkcji C.
        """
        assert type(wektor_1) == tuple, \
            'form of wektor_q is forbidden. wektor_1 should be touple'
        assert len(wektor_1) == 2,\
            'form of wektor_q is forbidden. wektor_1 should have two arguments'
        assert type(wektor_2) == tuple, \
            'form of wektor_q is forbidden. wektor_2 should be touple'
        assert len(wektor_2) == 2,\
            'form of wektor_q is forbidden. wektor_2 should have two arguments'
        assert type(znak) == str, \
            'znak is sign between two vector. Should be string'
        assert znak == '+' or znak == '-', \
            'only - and + are permitted'

        wekt_wypadkowy = self.norma_wektorow(wektor_1, wektor_2, znak)
        return cosh(wekt_wypadkowy * self.x) * exp(-wekt_wypadkowy * self.d / 2)

    def norma_wektorow(self, wektor_1, wektor_2, znak):
        """
        :param wektor_1: Pierwszy wektor do obliczenia różnicy.
        :param wektor_2: Drugi wektor do obliczenia różnicy.
        :param znak: Określa czy obliczana ma być różnica, czy suma wektorów.
        :return: W zlależności od znaku, zwraca normę z sumy, lub różnicy wektorów.
        """
        assert type(wektor_1) == tuple, \
            'form of wektor_q is forbidden. wektor_1 should be touple'
        assert len(wektor_1) == 2,\
            'form of wektor_q is forbidden. wektor_1 should have two arguments'
        assert type(wektor_2) == tuple, \
            'form of wektor_q is forbidden. wektor_2 should be touple'
        assert len(wektor_2) == 2,\
            'form of wektor_q is forbidden. wektor_2 should have two arguments'
        assert type(znak) == str, \
            'znak is sign between two vector. Should be string'
        assert znak == '+' or znak == '-', \
            'only - and + are permitted'

        wekt_wypadkowy = self.suma_roznica_wektorow(wektor_1, wektor_2, znak)
        return sqrt(wekt_wypadkowy[0] ** 2 + wekt_wypadkowy[1] ** 2)

    def delta_kroneckera(self):
        """
        Funkcja dodająca do macierzy pierwszy element z wyrażenia na M: 1 lub -1
        """
        for i in range(self.rozmiar_macierzy_blok, 2 * self.rozmiar_macierzy_blok):
            self.macierz_M[i - self.rozmiar_macierzy_blok][i] += 1
            self.macierz_M[i][i - self.rozmiar_macierzy_blok] -= 1

    def drugie_wyrazenie(self, wektor_1, wektor_2, wektor_q):
        """
        :param wektor_1: i-ty wektor
        :param wektor_2: j-ty wektor
        :param wektor_q: Blochowski wektor. Jest on "uciąglony". Jest on zmienną przy wyznaczaniu dyspersji.
        :return: Wynikiem jest drugie wyraz sumy.
        """
        assert type(wektor_1) == tuple, \
            'form of wektor_q is forbidden. wektor_1 should be touple'
        assert len(wektor_1) == 2,\
            'form of wektor_q is forbidden. wektor_1 should have two arguments'
        assert type(wektor_2) == tuple, \
            'form of wektor_q is forbidden. wektor_2 should be touple'
        assert len(wektor_2) == 2,\
            'form of wektor_q is forbidden. wektor_2 should have two arguments'
        assert type(wektor_q) == tuple, \
            'form of wektor_q is forbidden. wektor_q should be touple'
        assert len(wektor_q) == 2,\
            'form of wektor_q is forbidden. wektor_q should have two arguments'

        tmp1 = dot(self.suma_roznica_wektorow(wektor_q, wektor_2, "+"),
                   self.suma_roznica_wektorow(wektor_q, wektor_1, "+"))
        tmp2 = self.wektor_pozycji(wektor_1, wektor_2)
        return tmp1 * tmp2**2 / self.H0

    def trzecie_wyrazenie(self, wektor_1, wektor_2, wektor_q, typ_macierzy):
        """
        :param wektor_1: i-ty wektor.
        :param wektor_2: j-ty wektor.
        :param wektor_q: Blochowski wektor. Jest on "uciąglony". Jest on zmienną przy wyznaczaniu dyspersji.
        :param typ_macierzy: Określa, do której z macierzy blokowych odnosi się wyrażenie.
        :return: Wynikiem jest drugi wyraz sumy.
        """
        assert type(wektor_1) == tuple, \
            'form of wektor_q is forbidden. wektor_1 should be touple'
        assert len(wektor_1) == 2,\
            'form of wektor_q is forbidden. wektor_1 should have two arguments'
        assert type(wektor_2) == tuple, \
            'form of wektor_q is forbidden. wektor_2 should be touple'
        assert len(wektor_2) == 2,\
            'form of wektor_q is forbidden. wektor_2 should have two arguments'
        assert type(wektor_q) == tuple, \
            'form of wektor_q is forbidden. wektor_q should be touple'
        assert len(wektor_q) == 2,\
            'form of wektor_q is forbidden. wektor_q should have two arguments'
        assert typ_macierzy == 'xy' or typ_macierzy == 'yx',\
            'it is assumed that block matrixes are named xy or yx'

        tmp1 = self.norma_wektorow(wektor_q, wektor_2, '+')
        tmp2 = 1 - self.funkcja_c(wektor_q, wektor_2, "+")
        tmp3 = self.wspolczynnik(wektor_1, wektor_2)
        if typ_macierzy == 'xy':
            return (wektor_q[0] + wektor_2[0]) ** 2 / (self.H0 * tmp1 ** 2) * tmp2 * tmp3
        elif typ_macierzy == 'yx':
            return tmp2 * tmp3 / self.H0

    def czwarte_wyrazenie(self, wektor_1, wektor_2):
        """
        :param wektor_1: i-ty wektor.
        :param wektor_2: j-ty wektor.
        :return: Wynikiem jest czwarte wyrażenie w sumie na element macierzy M.
        """
        assert type(wektor_1) == tuple, \
            'form of wektor_q is forbidden. wektor_1 should be touple'
        assert len(wektor_1) == 2,\
            'form of wektor_q is forbidden. wektor_1 should have two arguments'
        assert type(wektor_2) == tuple, \
            'form of wektor_q is forbidden. wektor_2 should be touple'
        assert len(wektor_2) == 2,\
            'form of wektor_q is forbidden. wektor_2 should have two arguments'
        # TODO Poprawwić mianownik
        t = self.norma_wektorow(wektor_1, wektor_2, "-")
        if t == 0:
            tmp1 = 0
        else:
            tmp1 = (wektor_1[1] - wektor_2[1]) ** 2 / (self.H0 * t ** 2)
        tmp2 = self.wspolczynnik(wektor_1, wektor_2)
        tmp3 = (1 - self.funkcja_c(wektor_1, wektor_2, "-"))
        return tmp1 * tmp2 * tmp3

    def macierz_xy(self, wektor_1, wektor_2, wektor_q):
        assert type(wektor_1) == tuple, \
            'form of wektor_q is forbidden. wektor_1 should be touple'
        assert len(wektor_1) == 2,\
            'form of wektor_q is forbidden. wektor_1 should have two arguments'
        assert type(wektor_2) == tuple, \
            'form of wektor_q is forbidden. wektor_2 should be touple'
        assert len(wektor_2) == 2,\
            'form of wektor_q is forbidden. wektor_2 should have two arguments'
        assert type(wektor_q) == tuple, \
            'form of wektor_q is forbidden. wektor_q should be touple'
        assert len(wektor_q) == 2,\
            'form of wektor_q is forbidden. wektor_q should have two arguments'
        # TODO Dokończyć dokumentację
        return \
            self.drugie_wyrazenie(wektor_1, wektor_2, wektor_q) \
            + self.trzecie_wyrazenie(wektor_1, wektor_2, wektor_q, "xy") \
            - self.czwarte_wyrazenie(wektor_1, wektor_2)

    def macierz_yx(self, wektor_1, wektor_2, wektor_q):
        assert type(wektor_1) == tuple, \
            'form of wektor_q is forbidden. wektor_1 should be touple'
        assert len(wektor_1) == 2,\
            'form of wektor_q is forbidden. wektor_1 should have two arguments'
        assert type(wektor_2) == tuple, \
            'form of wektor_q is forbidden. wektor_2 should be touple'
        assert len(wektor_2) == 2,\
            'form of wektor_q is forbidden. wektor_2 should have two arguments'
        assert type(wektor_q) == tuple, \
            'form of wektor_q is forbidden. wektor_q should be touple'
        assert len(wektor_q) == 2,\
            'form of wektor_q is forbidden. wektor_q should have two arguments'

        # TODO Dokończyć dokumentację
        return \
            - self.drugie_wyrazenie(wektor_1, wektor_2, wektor_q) \
            - self.trzecie_wyrazenie(wektor_1, wektor_2, wektor_q, "yx") \
            + self.czwarte_wyrazenie(wektor_1, wektor_2)

    def lista_wektorow(self):
        # TODO Dokończyć dokumentację
        lista = WektorySieciOdwrotnej(self.a, self.a, self.rozmiar_macierzy_blok)
        return lista.lista_wektorow

    def wypelnienie_macierzy(self, wektor_q):
        assert type(wektor_q) == tuple, \
            'form of wektor_q is forbidden. wektor_q should be touple'
        assert len(wektor_q) == 2,\
            'form of wektor_q is forbidden. wektor_q should have two arguments'

        indeks = self.rozmiar_macierzy_blok
        lista_wektorow = self.lista_wektorow()
        assert len(lista_wektorow) == indeks ** 2, 'number of vector do not fit to matrix'
        self.delta_kroneckera()
        for i in range(indeks, 2 * indeks):
            for j in range(0, indeks):
                self.macierz_M[i][j] = \
                    self.macierz_xy(lista_wektorow[i - indeks], lista_wektorow[j], wektor_q)
                self.macierz_M[i - indeks][j + indeks] = \
                    self.macierz_yx(lista_wektorow[i], lista_wektorow[j - indeks], wektor_q)
        return self.macierz_M

    def wypisz_macierz(self):
        savetxt('macierz.txt', array(self.macierz_M))

# TODO rozdzielenie obliczeń dla wczytywania wektorów i współczynników z pliku oraz na klasę wykonującą te obliczenia analitycznie

        # q = MacierzDoZagadnienia(5)
        # print(q.czwarte_wyrazenie((-209439510.23931956, 418879020.4786391), (0.0, 209439510.23931956)))
        # print(q.wypelnienie_macierzy((3.807991095260355622e+07, 0)))
