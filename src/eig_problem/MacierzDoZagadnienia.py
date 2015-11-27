# -*- coding: utf-8 -*-


from numpy import sqrt, cosh, exp, dot
from numpy import zeros, array, savetxt

from src.eig_problem.FFTfromFile import FFTfromFile
from src.eig_problem.ParametryMaterialowe import ParametryMaterialowe
from src.eig_problem.WektorySieciOdwrotnej import WektorySieciOdwrotnej
from src.eig_problem.DFT import DFT

class MacierzDoZagadnienia(ParametryMaterialowe):
    """
    Klasa, w której tworzona jest macierz zagadnienia własnego. Pobiera ona współczynniki Fouriera z klasy FFTfromFile.
    Współczynniki są dla niej tworzone poprzez FFT z pliku graficznego.
    """

    def __init__(self, ilosc_wektorow, skad_wspolczynnik):
        ParametryMaterialowe.__init__(self, ilosc_wektorow)
        self.macierz_M = zeros((2 * ilosc_wektorow, 2 * ilosc_wektorow), dtype=complex)
        self.ilosc_wektorow = ilosc_wektorow
        if skad_wspolczynnik == 'FFT':
            self.slownik_wspolczynnik = FFTfromFile(ilosc_wektorow).fourier_coefficient()
            self.slownik_dlugosc_wymiany = FFTfromFile(ilosc_wektorow).exchange_length()
            self.lista_wektorow = FFTfromFile(ilosc_wektorow).vector_to_matrix()
        else:
            self.slownik_dlugosc_wymiany = DFT(self.ilosc_wektorow).slownik_wspolczynnikow()[1]
            self.slownik_wspolczynnik = DFT(self.ilosc_wektorow).slownik_wspolczynnikow()[0]
            self.lista_wektorow = WektorySieciOdwrotnej(self.a, self.b, ilosc_wektorow).lista_wektorow('min')

    def wspolczynnik(self, wektor_1, wektor_2):
        """
        Metoda wyliczająca współczynnik Fouriera. Jako argumenty podawne są dwa wektory i wyliczana jest
        z nich różnica.
        :type wektor_1: tuple
        :type wektor_2: tuple
        :param wektor_1: Pierwszy wektor do obliczenia różnicy.
        :param wektor_2: Drugi wektor do obliczenia różnicy.
        :return: współczynnik Fouriera dla różnicy wektorów sieci odwrotnej.
        """
        wekt_wypadkowy = self.suma_roznica_wektorow(wektor_1, wektor_2, '-')
        return self.slownik_wspolczynnik[wekt_wypadkowy]

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
        wekt_wypadkowy = self.suma_roznica_wektorow(wektor_1, wektor_2, '-')
        return self.slownik_dlugosc_wymiany[wekt_wypadkowy]

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
        assert len(wektor_1) == 2, \
            'form of wektor_q is forbidden. wektor_1 should have two arguments'
        assert len(wektor_2) == 2, \
            'form of wektor_q is forbidden. wektor_2 should have two arguments'
        assert type(znak) == str, \
            'znak is sign between two vector. Should be string'
        assert znak == '+' or znak == '-', \
            'only - and + are permitted'

        wekt_wypadkowy = self.norma_wektorow(wektor_1, wektor_2, znak)
        return cosh(wekt_wypadkowy * self.x) * exp(-wekt_wypadkowy * self.d / 2)

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
        assert len(wektor_1) == 2, \
            'form of wektor_q is forbidden. wektor_1 should have two arguments'
        assert len(wektor_2) == 2, \
            'form of wektor_q is forbidden. wektor_2 should have two arguments'
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
            self.macierz_M[i - self.ilosc_wektorow][i] += 1
            self.macierz_M[i][i - self.ilosc_wektorow] -= 1

    def drugie_wyrazenie(self, wektor_1, wektor_2, wektor_q):
        """
        Metoda obliczająca człon elmentu macierzowego związana z polem wymiany.
        :type wektor_1: tuple
        :type wektor_2: tuple
        :type wektor_q: tuple
        :param wektor_1: i-ty wektor
        :param wektor_2: j-ty wektor
        :param wektor_q: Blochowski wektor. Jest on "uciąglony". Jest on zmienną przy wyznaczaniu dyspersji.
        :return: Wynikiem jest drugie wyraz sumy.
        """
        assert len(wektor_1) == 2, \
            'form of wektor_q is forbidden. wektor_1 should have two arguments'
        assert len(wektor_2) == 2, \
            'form of wektor_q is forbidden. wektor_2 should have two arguments'
        assert len(wektor_q) == 2, \
            'form of wektor_q is forbidden. wektor_q should have two arguments'

        tmp1 = dot(self.suma_roznica_wektorow(wektor_q, wektor_2, "+"),
                   self.suma_roznica_wektorow(wektor_q, wektor_1, "+"))
        tmp2 = self.dlugosc_wymiany(wektor_1, wektor_2)
        return tmp1 * tmp2 / self.H0

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
        assert len(wektor_1) == 2, \
            'form of wektor_q is forbidden. wektor_1 should have two arguments'
        assert len(wektor_2) == 2, \
            'form of wektor_q is forbidden. wektor_2 should have two arguments'
        assert len(wektor_q) == 2, \
            'form of wektor_q is forbidden. wektor_q should have two arguments'
        assert typ_macierzy == 'xy' or typ_macierzy == 'yx', \
            'it is assumed that block matrixes are named xy or yx'
        tmp1 = self.norma_wektorow(wektor_q, wektor_2, '+')
        assert tmp1 != 0, 'probably insert forbidden q vector e.g. q = 0, q = 1'
        tmp2 = self.funkcja_c(wektor_q, wektor_2, "+")
        tmp3 = self.wspolczynnik(wektor_1, wektor_2)
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
        assert len(wektor_1) == 2, \
            'form of wektor_q is forbidden. wektor_1 should have two arguments'
        assert len(wektor_2) == 2, \
            'form of wektor_q is forbidden. wektor_2 should have two arguments'
        # TODO Poprawwić mianownik
        if self.norma_wektorow(wektor_1, wektor_2, "-") == 0:
            return 0
        else:
            return ((wektor_1[1] - wektor_2[1]) ** 2 / (self.H0 * self.norma_wektorow(wektor_1, wektor_2, "-") ** 2)) \
                   * self.wspolczynnik(wektor_1, wektor_2) \
                   * (1 - self.funkcja_c(wektor_1, wektor_2, "-"))

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
        # TODO Dokończyć dokumentację
        return \
            self.drugie_wyrazenie(wektor_1, wektor_2, wektor_q) \
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
            - self.drugie_wyrazenie(wektor_1, wektor_2, wektor_q) \
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
                self.macierz_M[i][j] += \
                    self.macierz_yx(lista_wektorow[i - indeks], lista_wektorow[j], wektor_q)
                self.macierz_M[i - indeks][j + indeks] += \
                    self.macierz_xy(lista_wektorow[i - indeks], lista_wektorow[j], wektor_q)
        return self.macierz_M

    def wypisz_macierz(self):
        """
        :return: Wypisuje tablice do pliku tekstowego.
         Ważne! Przed wypisaniem, należy wypełnić macierz_M metodą 'wypełnienie_macierzy'
        """
        savetxt('macierz.txt', array(self.macierz_M))

    if __name__ == 'main':
        pass
