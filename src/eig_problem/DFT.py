import numpy as np
from scipy import special

from src.eig_problem.ParametryMaterialowe import ParametryMaterialowe
from src.eig_problem.WektorySieciOdwrotnej import WektorySieciOdwrotnej


class DFT(ParametryMaterialowe):
    """
    Klasa obliczająca współczynniki Fouriera w sposób analityczny.
    """

    def __init__(self, ilosc_wektorow, typ_pole_wymiany):
        ParametryMaterialowe.__init__(self, ilosc_wektorow, typ_pole_wymiany)
        self.lista_wektorow = WektorySieciOdwrotnej(self.a, self.a, self.ilosc_wektorow).lista_wektorow('max')

    def wspolczynnik(self, wektor):
        """
        Metoda wyliczająca współczynnik Fouriera. Jako argumenty podawne są dwa wektory i wyliczana jest
        z nich różnica.
        :type wektor: tuple
        :param wektor: Wektor sieci odwrotnej, dla którego obliczany jest współczynnik Fouriera.
        :return: Współczynnik Fouriera dla zadanego wektora.
        """
        assert len(wektor) == 2, \
            'form of wektor_q is forbidden. wektor_1 should have two arguments'
        if wektor[0] == 0 and wektor[1] == 0:
            return (self.MoA - self.MoB) * np.pi * self.r ** 2 / (self.a ** 2) + self.MoB
        else:
            assert wektor[0] ** 2 + wektor[1] ** 2 != 0, 'division by 0'
            return 2 * (self.MoA - self.MoB) * np.pi * self.r ** 2 / self.a ** 2 * \
                   special.j1(np.sqrt(wektor[0] ** 2 + wektor[1] ** 2) * self.r) / \
                   (np.sqrt(wektor[0] ** 2 + wektor[1] ** 2) * self.r)

    def dlugosc_wymiany(self, wektor):
        """
        Metoda obliczająca długość wymiany, dla dwóch zadanych wektorów.
        z nich różnica.
        :type wektor: tuple
        :param wektor: Wektor sieci odwrotnej, dla którego obliczany jest współczynnik Fouriera.
        :return: Długość wymiany dla zadanego wektora.
        """
        assert len(wektor) == 2, \
            'form of wektor_q is forbidden. wektor_1 should have two arguments'

        if wektor[0] == 0 and wektor[1] == 0:
            return (self.lA - self.lB) * np.pi * self.r ** 2 / (self.a ** 2) + self.lB
        else:
            return 2 * (self.lA - self.lB) * np.pi * self.r ** 2 / (self.a ** 2) * \
                   special.j1(np.sqrt(wektor[0] ** 2 + wektor[1] ** 2) * self.r) \
                   / (np.sqrt(wektor[0] ** 2 + wektor[1] ** 2) * self.r)

    def slownik_wspolczynnikow(self):
        """
        Metoda, której zadaniem jest odpowiednie zestawienie współczynników Fouriera oraz długości wymiany
        wraz z wektorami sieci odwrotnej.
        :return: Zwraca dwa słowniki. Odpowiednio wektory sieci odwrotnej ze współczynnikami oraz wektory sieci
        odwrotnej z długościami wymiany.
        """
        lista_wektorow = self.lista_wektorow
        dlugosc_wymiany = np.zeros(len(lista_wektorow), dtype=complex)
        wspolczynnik = np.zeros(len(lista_wektorow), dtype=complex)
        for i, j in list(enumerate(lista_wektorow)):
            dlugosc_wymiany[int(i)] = self.dlugosc_wymiany(j)
            wspolczynnik[int(i)] = self.wspolczynnik(j)
        return dict(list(zip(lista_wektorow, wspolczynnik))), dict(list(zip(lista_wektorow, dlugosc_wymiany)))
