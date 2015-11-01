from numpy import zeros
from numpy import pi, sqrt, cosh, linalg, exp
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

    def __init__(self, rozmiar_macierzy_blok):
        self.macierz_M = zeros((2 * rozmiar_macierzy_blok, 2 * rozmiar_macierzy_blok))
        self.rozmiar_macierzy_blok = rozmiar_macierzy_blok

    def delta_kroneckera(self):
        """
        Funkcja dodająca do macierzy pierwszy element z wyrażenia na M: 1 lub -1
        """
        for i in range(self.rozmiar_macierzy_blok, 2 * self.rozmiar_macierzy_blok):
            self.macierz_M[i - self.rozmiar_macierzy_blok][i] += 1
            self.macierz_M[i][i - self.rozmiar_macierzy_blok] -= 1

    def wspolczynnik(self, wektor_g1, wektor_g2):
        """
        Metoda wyliczająca współczynnik Fouriera. Jako argumenty podawne są dwa wektory i wyliczana jest
        z nich różnica.
        :param wektor_g1: Pierwszy wektor do obliczenia różnicy.
        :param wektor_g2: Drugi wektor do obliczenia różnicy.
        :return: współczynnik Fouriera.
        """
        zipped = list(zip(wektor_g1, wektor_g2))
        wekt_wypadkowy = [k[0] - k[1] for k in zipped]

        if wekt_wypadkowy[0] == 0 and wekt_wypadkowy[1] == 0:
            return (self.MoA - self.MoB) * pi * self.r ** 2 / self.s + self.MoB
        else:
            return 2 * (self.MoA - self.MoB) * pi * self.r ** 2 / self.s * \
                   scipy.special.j1(sqrt(wekt_wypadkowy[0] ** 2 + wekt_wypadkowy[1] ** 2) * self.r) \
                   / (sqrt(wekt_wypadkowy[0] ** 2 + wekt_wypadkowy[1] ** 2 + (10 ** -10)) * self.r)

    @staticmethod
    def sprawdz_znak(wektor_1, wektor_2, znak):
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
        wekt_wypadkowy = self.sprawdz_znak(wektor_1, wektor_2, znak)
        return cosh(linalg.norm(wekt_wypadkowy) * x) * exp(-linalg.norm(wekt_wypadkowy) * d / 2)

    def norma_wektorow(self, wektor_1, wektor_2, znak):
        """
        :param wektor_1: Pierwszy wektor do obliczenia różnicy.
        :param wektor_2: Drugi wektor do obliczenia różnicy.
        :param znak: Określa czy obliczana ma być różnica, czy suma wektorów.
        :return: W zlależności od znaku, zwraca normę z sumy, lub różnicy wektorów.
        """
        wekt_wypadkowy = self.sprawdz_znak(wektor_1, wektor_2, znak)
        return linalg.norm(wekt_wypadkowy)

    def czwarte_wyrazenie(self):

        # mxy = zeros((self.rozmiar_macierzy_blok, self.rozmiar_macierzy_blok))
        myx = zeros((self.rozmiar_macierzy_blok, self.rozmiar_macierzy_blok))

        return myx

    def wypisz_macierz(self):
        print(self.macierz_M)


b = MacierzDoZagadnienia(5)
# print(b.wspolczynnik((1, 1), (0, 0)))
# print(b.wektory(1))
print(b.funkcja_c((1, 0), (3, 0), 3, 5, "-"))
# print(b.norma_wektorow(3,4, "+"))
b.delta_kroneckera()
b.wypisz_macierz()
