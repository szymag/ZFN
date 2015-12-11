import cProfile

from numpy import linspace, pi, savetxt
from scipy import linalg

from src.eig_problem.MacierzDoZagadnienia import MacierzDoZagadnienia
from src.eig_problem.ParametryMaterialowe import ParametryMaterialowe


def do_cprofile(func):
    def profiled_func(*args, **kwargs):
        profile = cProfile.Profile()
        try:
            profile.enable()
            result = func(*args, **kwargs)
            profile.disable()
            return result
        finally:
            profile.print_stats()

    return profiled_func


# @do_cprofile
class ZagadnienieWlasne(ParametryMaterialowe):
    """
    Klasa, w której zdefiniowane jest uogólnione zagadnienie własne. Dziedziczy ona po 'Parametry materiałowe, gdyż
    potrzebne są w niej infoemacje o strukturze kryształu magnonicznego.
    """
    def __init__(self, ilosc_wektorow, ilosc_wektorow_q, skad_wspolczynnik, typ_pole_wymiany):
        """
        :param ilosc_wektorow: Ile wektorów wchodzi do zagadnienia własnego. Determinuje to wielkość macierzy blokowych.
        :param ilosc_wektorow_q: Odpowiada za gęstość siatki, na wykresie dyspersji.
        """
        ParametryMaterialowe.__init__(self, ilosc_wektorow, typ_pole_wymiany)
        self.lista_wektorow_q = [((2 * pi * k / self.a), pi / self.a) for k in linspace(0.01, 0.5, ilosc_wektorow_q)]
        self.skad_wspolczynnik = skad_wspolczynnik
        self.typ_pola_wymiany = typ_pole_wymiany

    # @do_cprofile
    def utworz_macierz_m(self, wektor_q):
        """
        :rtype: object
        :type wektor_q tuple
        :param wektor_q: Blochowski wektor. Jest on "uciąglony". Jest on zmienną przy wyznaczaniu dyspersji.
        :return: Zwraca tablicę do zagadnienia wlasnego, która została utworzona w klasie "MacierzDoZagadnienia".
        """
        assert type(wektor_q) == tuple, \
            'form of wektor_q is forbidden. wektor_q should be touple'
        assert len(wektor_q) == 2,\
            'form of wektor_q is forbidden. wektor_q should have two arguments'
        macierz = MacierzDoZagadnienia(self.ilosc_wektorow, self.skad_wspolczynnik,
                                       self.typ_pola_wymiany).wypelnienie_macierzy(wektor_q)
        return macierz

    def zagadnienie_wlasne(self, wektor_q):
        """
        Metoda, która wywołuje algorytm rozwiązywania zagadnienia własnego. Poprzez metodę 'utworz_macierz_M' tworzy
        sobie tablicę, dla której następnie oblicza wartości i wektory własne. Jest ona także przystosowana dla
        uogólnionego zagadnienia własnego.
        :type wektor_q tuple
        :param wektor_q: Blochowski wektor. Jest on "uciąglony". Jest on zmienną przy wyznaczaniu dyspersji.
        :return: Wartości własne. Wektory własne są obecnie wyłączone.
        """
        # m_jednostkowa = identity(2 * self.ilosc_wektorow)
        assert type(wektor_q) == tuple, \
            'form of wektor_q is forbidden. wektor_q should be touple'
        assert len(wektor_q) == 2,\
            'form of wektor_q is forbidden. wektor_q should have two arguments'
        # TODO: słówko yield!
        macierz_m = self.utworz_macierz_m(wektor_q)
        return linalg.eig(macierz_m, right=False)

    def czestosci_wlasne(self, wektor_q):
        """
        Metoda, w której przelicza się wartości własne w metodzie 'zagadnienie wlasne' na czestosci wlasne' wzbudzen
        w krysztale magnonicznym.
        :type wektor_q: tuple
        :param wektor_q:  Blochowski wektor. Jest on "uciąglony". Jest on zmienną przy wyznaczaniu dyspersji.
        :return: Częstości własne, obliczane według wzoru: w = -a.imag * gamma * mu0*H0, gdzie w - cżęstość własne,
        a.imag część urojona wartości własnej. Ponieważ wartości własne są przmnażane przez 1/i, wystarczy wziąć część
        zespoloną, by uzysać częstotliwość rzeczywistą. Część zespolona w częstotliwośći odpowiada za tłumienie
        w układzie, a takiego się nie zakłada, więc powinny one być stosunkowo niewielkie - pomijalne.
        """
        assert type(wektor_q) == tuple, \
            'form of wektor_q is forbidden. wektor_q should be touple'
        assert len(wektor_q) == 2,\
            'form of wektor_q is forbidden. wektor_q should have two arguments'
        wartosci_wlasne = self.zagadnienie_wlasne(wektor_q)
        czestosci_wlasne = [i.imag * self.gamma * self.mu0H0 / 2 / pi for i in wartosci_wlasne]
        czestosci_wlasne = [i for i in czestosci_wlasne if i > 0]
        czestosci_wlasne = list(sorted(czestosci_wlasne))
        return czestosci_wlasne[0:10]

    def wypisz_do_pliku(self):
        """
        Metoda, która wypisuje w kolumnach wyniki. W pierwszej jest wektor blochowski, a w dalszych odpowiadające mu
        częstotliwości własne.
        :return: Plik tekstowy.
        """
        plik = []
        for k in self.lista_wektorow_q:
            tmp = [k[0]]
            tmp.extend(self.czestosci_wlasne(k))
            plik.append(tmp)
        savetxt('1.txt', plik)


def start():
    return ZagadnienieWlasne(225, 20, 'FFT', 'II').wypisz_do_pliku()


start()
