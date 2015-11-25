import numpy as np
from scipy import special

from src.eig_problem.ParametryMaterialowe import ParametryMaterialowe
from src.eig_problem.WektorySieciOdwrotnej import WektorySieciOdwrotnej


class DFT(ParametryMaterialowe):
    def __init__(self, ilosc_wektorow):
        ParametryMaterialowe.__init__(self, ilosc_wektorow)

    @staticmethod
    def suma_roznica_wektorow(wektor_1, wektor_2, znak):
        """
        Metoda, która w zależności od znaku oblicza sumę, bądż różnicę wektorów.
        :param wektor_1: Pierwszy wektor do obliczenia różnicy.
        :param wektor_2: Drugi wektor do obliczenia różnicy.
        :param znak: Określa czy obliczana ma być różnica, czy suma wektorów.
        :return: suma lub rożnica wektorów
        """
        assert len(wektor_1) == 2, \
            'form of wektor_q is forbidden. wektor_1 should have two arguments'
        assert len(wektor_2) == 2, \
            'form of wektor_q is forbidden. wektor_2 should have two arguments'
        assert type(znak) == str, \
            'znak is sign between two vector. Should be string'
        assert znak == '+' or znak == '-', \
            'only - and + are permitted'

        if znak == "-":
            return tuple([k[0] - k[1] for k in zip(wektor_1, wektor_2)])
        elif znak == "+":
            return tuple([k[0] + k[1] for k in zip(wektor_1, wektor_2)])

    def lista_wektorow(self):
        # TODO Dokończyć dokumentację
        indeks = self.ilosc_wektorow
        lista = WektorySieciOdwrotnej(self.a, self.a, indeks)
        return lista.lista_wektorow()

    def wspolczynnik(self, wektor_1, wektor_2):
        """
        Metoda wyliczająca współczynnik Fouriera. Jako argumenty podawne są dwa wektory i wyliczana jest
        z nich różnica.
        :param wektor_1: Pierwszy wektor do obliczenia różnicy.
        :param wektor_2: Drugi wektor do obliczenia różnicy.
        :return: współczynnik Fouriera dla różnicy wektorów sieci odwrotnej.
        """
        assert len(wektor_1) == 2, \
            'form of wektor_q is forbidden. wektor_1 should have two arguments'
        assert len(wektor_2) == 2, \
            'form of wektor_q is forbidden. wektor_2 should have two arguments'
        wekt_wypadkowy = self.suma_roznica_wektorow(wektor_1, wektor_2, '-')
        if wekt_wypadkowy[0] == 0 and wekt_wypadkowy[1] == 0:
            return (self.MoCo - self.MoPy) * pi * self.r ** 2 / (self.a ** 2) + self.MoPy
        else:
            assert wekt_wypadkowy[0] ** 2 + wekt_wypadkowy[1] ** 2 != 0, 'division by 0'
            return 2 * (self.MoCo - self.MoPy) * pi * self.r ** 2 / self.a ** 2 * \
                   special.j1(np.sqrt(wekt_wypadkowy[0] ** 2 + wekt_wypadkowy[1] ** 2) * self.r) / \
                   (np.sqrt(wekt_wypadkowy[0] ** 2 + wekt_wypadkowy[1] ** 2) * self.r)

    def dlugosc_wymiany(self, wektor_1, wektor_2):
        """
        Metoda obliczająca długość wymiany, dla dwóch zadanych wektorów.
        z nich różnica.
        :param wektor_1: Pierwszy wektor do obliczenia różnicy.
        :param wektor_2: Drugi wektor do obliczenia różnicy.
        :return: Długość wymiany w postaci odpowiadającego różnicy wektorów współczynnika Fouriera.
        """
        assert len(wektor_1) == 2, \
            'form of wektor_q is forbidden. wektor_1 should have two arguments'
        assert len(wektor_2) == 2, \
            'form of wektor_q is forbidden. wektor_2 should have two arguments'

        wekt_wypadkowy = self.suma_roznica_wektorow(wektor_1, wektor_2, '-')

        if wekt_wypadkowy[0] == 0 and wekt_wypadkowy[1] == 0:
            return (self.lCo - self.lPy) * pi * self.r ** 2 / (self.a ** 2) + self.lPy
        else:
            return 2 * (self.lCo - self.lPy) * pi * self.r ** 2 / (self.a ** 2) * \
                   special.j1(np.sqrt(wekt_wypadkowy[0] ** 2 + wekt_wypadkowy[1] ** 2) * self.r) \
                   / (np.sqrt(wekt_wypadkowy[0] ** 2 + wekt_wypadkowy[1] ** 2) * self.r)
