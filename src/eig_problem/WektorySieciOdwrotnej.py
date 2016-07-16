from math import sqrt

import numpy as np


class WektorySieciOdwrotnej:
    """
    Klasa, której zadaniem jest wygenerowanie wektorow sieci odwrotnej, służących dalej do definiowania
    macierzy zagadnienia własnego. W tym przypadku wykorzystwyane są wzory analityczne na współczynniki Fouriera.
    """

    def __init__(self, dlugosc_skladowej_y, dlugosc_skladowej_z, rozmiar_macierzy_blok):
        self.dlugosc_skladowej_y = dlugosc_skladowej_y
        self.dlugosc_skladowej_z = dlugosc_skladowej_z
        self.ilosc_wektorow = rozmiar_macierzy_blok

    def wspolrzedna_wektora(self, k, typ):
        """
        Metoda zwracająca listę współrzędnych wektora sieci odwrotnej zdefiniowanego następująco: 2*pi*n/a, gdzie a jest
        stałą sieciową, a n liczbą całkowitą numerującą wektor w danym zakresie(-nmax, nmax). Dodatkowo, na metodę
        narzucone są warunki dotyczące liczby wektorów. Musi ich być nieprzyście oraz muszą być pierwiastkiem jakieś
        liczby całkowitej.

        :param k: Określa, którą współrzędną metoda ma obliczyć. Dla k=1 igrekową, dla k=2 zetową.
        :param typ: Określa, z jakiego zakresu mają być wybierane wektory sieci odwrotnej. Wartość min oznacza, że
        z podstawowego tzn. dla tych wektorów generuje się zagadnienie własne. Wartość max określa zakres dwa razy
        większy, by obliczone wartości pewnych parametrów także zostały obliczone. Np. różnica dwóch wektórów
        maksymalnych z zakresu min daje wartość 2*mi
        :return: Zwraca odpowiednią listę współrzędnych wektorów sieci odwrotnej.
        """
        assert k == 1 or k == 2, 'k should bo 1 or 2'
        assert sqrt(self.ilosc_wektorow) == int(sqrt(self.ilosc_wektorow)) and self.ilosc_wektorow % 2 != 0, \
            'size of bolck matrix shuld gave natural number of sqrt and be odd'
        assert typ == 'max' or typ == 'min', 'typ should by max or min'
        if typ == 'max':
            indeks = int(2 * sqrt(self.ilosc_wektorow) / 2) - 1
        else:
            indeks = int(sqrt(self.ilosc_wektorow + 1) / 2)
        lista = np.array(range(-indeks, indeks + 1))
        if k == 1:
            lista = [int(2 * 3 * i) for i in lista]
        elif k == 2:
            lista = [int(2 * 3 * i) for i in lista]
        return lista

    def lista_wektorow(self, typ):
        """
        Metoda obliczająca iloczyn kartezjański dla dwóch wygenerowanych list metodą 'wspolrzedna_wektora'.
        :param typ: Określa, z jakiego zakresu mają być wybierane wektory sieci odwrotnej. Wartość min oznacza, że
        z podstawowego tzn. dla tych wektorów generuje się zagadnienie własne. Wartość max określa zakres dwa razy
        większy, by obliczone wartości pewnych parametrów także zostały obliczone. Np. różnica dwóch wektórów
        maksymalnych z zakresu min daje wartość 2*min.
        :return: Lista wektorów sieci odwrotnej.
        """
        lista = []
        lista1 = self.wspolrzedna_wektora(1, typ)
        lista2 = self.wspolrzedna_wektora(2, typ)
        for i in lista1:
            for j in lista2:
                lista.append((i, j))
        # lista = list(product(*(self.wspolrzedna_wektora(1, typ), self.wspolrzedna_wektora(2, typ))))
        return lista

    def wspolrzedna_wektora1d(self, typ):
        assert typ == 'max' or typ == 'min', 'typ should by max or min'
        if typ == 'min':
            return np.arange(-(self.ilosc_wektorow - 1) // 2, (self.ilosc_wektorow - 1) // 2 + 1)
        elif typ == 'max':
            return np.arange(- (self.ilosc_wektorow - 1), (self.ilosc_wektorow - 1) + 1)

    def lista_wektorow1d(self, typ):
        return self.wspolrzedna_wektora1d(typ)

