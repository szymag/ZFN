import math

import numpy as np


class WektorySieci(object):
    """
    Klasa, kótrej zadaniem jest wyznaczenie wektorów sieci odwrotnej dla zadanych wektorów sieci rzeczywistej
    """
    a_z = [0, 0, 1]

    def __init__(self, dlugosc_a1, dlugosc_a2, kat, zakres_1, zakres_2):
        """

        :param dlugosc_a1: długość wektora a1 sieci rzeczywistej.
        :param dlugosc_a2: długość wektora a2 sieci rzeczywistej
        :param kat: określa, kąt pomiędzy dwoma wektorami sieci rzeczywistej
        :param zakres_1: określa, ile wektorów sieci odwrotnej b1 zostanie wygenerowanych
        :param zakres_2: określa, ile wektorów sieci odwrotnej b2 zostanie wygenerowanych
        """
        self.dlugosc_a1 = dlugosc_a1
        self.dlugosc_a2 = dlugosc_a2
        self.kat = kat
        self.zakres_1 = list(range(-zakres_1, zakres_1 + 1))
        self.zakres_2 = list(range(-zakres_2, zakres_2 + 1))

    def wektor_a1(self):
        """
        :return: funkcja zwraca listę z współrzędnymi wektora a1
        """
        return [self.dlugosc_a1, 0, 0]

    def wektor_a2(self):
        """
        :return: funkcja zwraca listę z współrzędnymi wektora a2
        """
        return [self.dlugosc_a2 * math.cos(math.radians(self.kat)),
                self.dlugosc_a2 * math.sin(math.radians(self.kat)),
                0]

    def wektor_b(self, k):
        """
        :param k: określa, który wektor, b1 czy b2 zostanie obliczony.
        :return: funkcja zwraca listę z współrzędnymi wektora bk gdzie k przyjmuje wartości 1 lub 2
        """
        a1 = self.wektor_a1()
        a2 = self.wektor_a2()
        n = self.a_z
        temp = 2 * np.pi
        if k == 1:
            temp *= np.cross(a2, n)
        elif k == 2:
            temp *= np.cross(n, a1)
        return temp / np.dot(a1, np.cross(a2, n))

    def lista_wektorow_b1(self):
        return [self.wektor_b(1) * k for k in self.zakres_1]

    def lista_wektorow_b2(self):
        return [self.wektor_b(2) * k for k in self.zakres_2]

    def wektor_wypadkowy(self):
        return [a + b for a, b in zip(self.wektor_b(1), self.wektor_b(2))]

    def wektor_g(self, k):
        return [self.wektor_b(k)[k - 1] * i for i in self.zakres_1]
