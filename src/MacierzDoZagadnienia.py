from numpy import zeros
from numpy import pi, sqrt, cosh, linalg, exp, linspace, dot
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

    def drugie_wyrazenie(self, wektor_1, wektor_2, wektor_q, typ=None):
        tmp1 = dot(self.suma_roznica_wektorow(wektor_q, wektor_2, "+"),
                   self.suma_roznica_wektorow(wektor_q, wektor_1, "+"))
        tmp2 = self.wspolczynnik(wektor_1, wektor_2)
        return tmp1 * tmp2 / self.H0

    def trzecie_wyrazenie(self):
        pass

    def czwarte_wyrazenie(self, wektor_1, wektor_2):
        """
        :param wektor_1: i-ty wektor
        :param wektor_2: j-ty wektor
        :return: Wynikiem jest czwarte wyrażenie w sumie na element macierzy M
        """
        tmp1 = (wektor_1[1] - wektor_2[1]) ** 2 / (self.H0 * self.norma_wektorow(wektor_1, wektor_2, "-"))
        tmp2 = self.wspolczynnik(wektor_1, wektor_2)
        tmp3 = (1 - self.funkcja_c(wektor_1, wektor_2, self.x, self.d, "-"))
        return tmp1 * tmp2 * tmp3

    def wektor(self):
        return self.lista_wektorow_q

    def wypisz_macierz(self):
        print(self.macierz_M)


b = MacierzDoZagadnienia(5, 3)
print(b.wspolczynnik((1, 1), (0, 0)))
print(b.funkcja_c((1, 0), (3, 0), 3, 5, "-"))
print(b.norma_wektorow((3, 5), (4, 4), "+"))
b.delta_kroneckera()
b.wypisz_macierz()
print(b.wektor())
