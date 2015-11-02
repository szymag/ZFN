from numpy import zeros
from numpy import pi, sqrt, cosh, linalg, exp, linspace, dot
from src.WektorySieci import WektorySieci
import scipy.special


class MacierzDoZagadnienia:
    mu0H0 = 0.001
    H0 = mu0H0 / (4 * pi * 10 ** -7)
    d = 6
    a = 10
    s = a ** 2
    MoA = 10
    MoB = 4
    r = 3
    x = 1

    def __init__(self, rozmiar_macierzy_blok, gestosc_wektorow_q):
        self.macierz_M = zeros((2 * rozmiar_macierzy_blok, 2 * rozmiar_macierzy_blok))
        self.rozmiar_macierzy_blok = rozmiar_macierzy_blok
        self.lista_wektorow_q = [((2 * pi * k / self.a), 0) for k in linspace(0, 1, gestosc_wektorow_q)]

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
            return (self.MoA - self.MoB) * pi * self.r ** 2 / self.s + self.MoB
        else:
            return 2 * (self.MoA - self.MoB) * pi * self.r ** 2 / self.s * \
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
        :return: Wynikiem jest drugie wyraz sumy.
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
        tmp1 = (wektor_1[1] - wektor_2[1]) ** 2 / (self.H0 * self.norma_wektorow(wektor_1, wektor_2, "-"))
        tmp2 = self.wspolczynnik(wektor_1, wektor_2)
        tmp3 = (1 - self.funkcja_c(wektor_1, wektor_2, self.x, self.d, "-"))
        return tmp1 * tmp2 * tmp3

    def macierz_xy(self, wektor_1, wektor_2, wektor_q):
        return \
            self.drugie_wyrazenie(wektor_1, wektor_2, wektor_q) \
            + self.trzecie_wyrazenie(wektor_1, wektor_2, wektor_q, "xy") \
            - self.czwarte_wyrazenie(wektor_1, wektor_2)

    def macierz_yx(self, wektor_1, wektor_2, wektor_q):
        return \
            - self.drugie_wyrazenie(wektor_1, wektor_2, wektor_q) \
            - self.trzecie_wyrazenie(wektor_1, wektor_2, wektor_q, "yx") \
            + self.czwarte_wyrazenie(wektor_1, wektor_2)

    def lista_wektorow(self):
        # TODO PRzerobić klasę WektorySieci, tak by mieć siatkę wektorów
        zakres = int((self.rozmiar_macierzy_blok - 1) / 2)
        lista = WektorySieci(self.a, self.a, 90, zakres, zakres).wektor_g()
        return lista

    def wypelnienie_macierzy(self, wektor_q):
        # TODO Dokończyć budowanie metody
        indeks = self.rozmiar_macierzy_blok
        self.delta_kroneckera()
        for i in range(indeks, 2 * indeks):
            for j in range(0, indeks):
                # self.macierz_M[i][j] = self.macierz_xy()
                pass

    def wektor(self):
        return self.lista_wektorow_q

    def wypisz_macierz(self):
        print(self.macierz_M)


b = MacierzDoZagadnienia(10, 3)
# b.wypelnienie_macierzy()
b.wypisz_macierz()

print(b.lista_wektorow())
