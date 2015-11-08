import scipy.special
from numpy import pi, sqrt, cosh, linalg, exp, dot
from numpy import zeros

from src.ParametryMaterialowe import ParametryMaterialowe
from src.WektorySieciOdwrotnej import WektorySieciOdwrotnej


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
        zipped = list(zip(wektor_1, wektor_2))
        wekt_wypadkowy = [k[0] + k[1] for k in zipped]

        if wekt_wypadkowy[0] == 0 and wekt_wypadkowy[1] == 0:
            return (self.MoA - self.MoB) * pi * self.r ** 2 / (self.a) ** 2 + self.MoB
        else:
            return 2 * (self.MoA - self.MoB) * pi * self.r ** 2 / (self.a) ** 2 * \
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
        zipped = list(zip(wektor_1, wektor_2))
        if znak == "-":
            return [k[0] - k[1] for k in zipped]
        elif znak == "+":
            return [k[0] + k[1] for k in zipped]

    def funkcja_c(self, wektor_1, wektor_2, x, d, znak):
        """
        Metoda obliczająca wartość funkcji C dla rożnicy wektorów
        :param wektor_1: Pierwszy wektor do obliczenia różnicy.
        :param wektor_2: Drugi wektor do obliczenia różnicy.
        :param x: Współrzędna x (wzdłuż grubości warstwy).
        :param d: grubość warstwy.
        :param znak: określa czy obliczana ma być różnica, czy suma wektorów.
        :return: wartość funkcji C.
        """
        wekt_wypadkowy = self.suma_roznica_wektorow(wektor_1, wektor_2, znak)
        return cosh(linalg.norm(wekt_wypadkowy) * x) * exp(-linalg.norm(wekt_wypadkowy) * d / 2)

    def norma_wektorow(self, wektor_1, wektor_2, znak):
        """
        :param wektor_1: Pierwszy wektor do obliczenia różnicy.
        :param wektor_2: Drugi wektor do obliczenia różnicy.
        :param znak: Określa czy obliczana ma być różnica, czy suma wektorów.
        :return: W zlależności od znaku, zwraca normę z sumy, lub różnicy wektorów.
        """
        wekt_wypadkowy = self.suma_roznica_wektorow(wektor_1, wektor_2, znak)
        return linalg.norm(wekt_wypadkowy)

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
        tmp1 = dot(self.suma_roznica_wektorow(wektor_q, wektor_2, "+"),
                   self.suma_roznica_wektorow(wektor_q, wektor_1, "+"))
        tmp2 = self.wspolczynnik(wektor_1, wektor_2)
        return tmp1 * tmp2 / self.H0

    def trzecie_wyrazenie(self, wektor_1, wektor_2, wektor_q, typ_macierzy):
        """
        :param wektor_1: i-ty wektor.
        :param wektor_2: j-ty wektor.
        :param wektor_q: Blochowski wektor. Jest on "uciąglony". Jest on zmienną przy wyznaczaniu dyspersji.
        :param typ_macierzy: Określa, do której z macierzy blokowych odnosi się wyrażenie.
        :return: Wynikiem jest drugi wyraz sumy.
        """
        tmp1 = self.norma_wektorow(wektor_q, wektor_2, "+") ** 2
        tmp2 = self.funkcja_c(wektor_q, wektor_2, self.x, self.d, "+")
        tmp3 = self.wspolczynnik(wektor_1, wektor_2)
        if typ_macierzy == "xy":
            return (wektor_q[0] + wektor_2[0]) ** 2 / (self.H0 * tmp1) * (1 - tmp2) * tmp3
        elif typ_macierzy == "yx":
            return tmp2 * tmp3 / self.H0

    def czwarte_wyrazenie(self, wektor_1, wektor_2):
        """
        :param wektor_1: i-ty wektor.
        :param wektor_2: j-ty wektor.
        :return: Wynikiem jest czwarte wyrażenie w sumie na element macierzy M.
        """
        tmp1 = (wektor_1[1] - wektor_2[1]) ** 2 / (10 ** -7 + self.H0 * self.norma_wektorow(wektor_1, wektor_2, "-"))
        tmp2 = self.wspolczynnik(wektor_1, wektor_2)
        tmp3 = (1 - self.funkcja_c(wektor_1, wektor_2, self.x, self.d, "-"))
        return tmp1 * tmp2 * tmp3

    def macierz_xy(self, wektor_1, wektor_2, wektor_q):
        # TODO Dokończyć dokumentację
        return \
            self.drugie_wyrazenie(wektor_1, wektor_2, wektor_q) \
            + self.trzecie_wyrazenie(wektor_1, wektor_2, wektor_q, "xy") \
            - self.czwarte_wyrazenie(wektor_1, wektor_2)

    def macierz_yx(self, wektor_1, wektor_2, wektor_q):
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
        # TODO Dokończyć budowanie metody
        indeks = self.rozmiar_macierzy_blok
        lista_wektorow = self.lista_wektorow()
        self.delta_kroneckera()
        for i in range(indeks, 2 * indeks):
            for j in range(0, indeks):
                self.macierz_M[i][j] = \
                    self.macierz_xy(lista_wektorow[i - indeks], lista_wektorow[j], wektor_q)
                self.macierz_M[i - indeks][j + indeks] = \
                    self.macierz_yx(lista_wektorow[i], lista_wektorow[j - indeks], wektor_q)
        return self.macierz_M

    def wypisz_macierz(self):
        print(self.macierz_M)
