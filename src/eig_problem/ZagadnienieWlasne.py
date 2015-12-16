from src.eig_problem.cProfiler import do_cprofile
from numpy import linspace, pi, savetxt, array
from scipy import linalg
from src.eig_problem.MacierzDoZagadnienia import MacierzDoZagadnienia
from src.eig_problem.ParametryMaterialowe import ParametryMaterialowe
from operator import itemgetter

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
        self.lista_wektorow_q = [((2 * pi * k / self.a), 0.0) for k in linspace(0.01, 0.99, ilosc_wektorow_q)]
        self.skad_wspolczynnik = skad_wspolczynnik
        self.typ_pola_wymiany = typ_pole_wymiany

    #@do_cprofile
    def zagadnienie_wlasne(self, wektor_q, param):
        """
        Metoda, która wywołuje algorytm rozwiązywania zagadnienia własnego. Tworzy sobie tablicę,
        dla której następnie oblicza wartości i wektory własne. Jest ona także przystosowana dla
        uogólnionego zagadnienia własnego.
        :type wektor_q tuple
        :param wektor_q: Blochowski wektor. Jest on "uciąglony". Jest on zmienną przy wyznaczaniu dyspersji.
        :return: Wartości własne. Wektory własne są obecnie wyłączone.
        """
        # TODO: słówko yield!
        macierz_m = MacierzDoZagadnienia(self.ilosc_wektorow, self.skad_wspolczynnik,
                                         self.typ_pola_wymiany).wypelnienie_macierzy(wektor_q)
        return linalg.eig(macierz_m, right=param)  # trzeba pamiętać o włączeniu/wyłączeniu generowania wektorów

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
        assert len(wektor_q) == 2, \
            'form of wektor_q is forbidden. wektor_q should have two arguments'
        wartosci_wlasne = self.zagadnienie_wlasne(wektor_q, param=False)
        # wartosci_wlasne, wektory_wlasny = self.zagadnienie_wlasne(wektor_q)
        # savetxt(str(wektor_q), wektory_wlasny.view(float))
        czestosci_wlasne = [i.imag * self.gamma * self.mu0H0 / 2.0 / pi for i in wartosci_wlasne if i.imag > 0]

        return list(sorted(czestosci_wlasne))

    def wypisz_czestosci_do_pliku(self):
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

    def wektory_wlasne(self):
        assert len(self.lista_wektorow_q) == 1
        wartosci_wlasne, wektory_wlasne = self.zagadnienie_wlasne(self.lista_wektorow_q[0], param=True)

        czestosci_wlasne = [i.imag * self.gamma * self.mu0H0 / 2.0 / pi for i in wartosci_wlasne if i.imag > 0]
        czestosci_wlasne = enumerate(czestosci_wlasne)
        czestosci_wlasne = list(sorted(czestosci_wlasne, key=itemgetter(1)))
        wektory_wlasne = list(enumerate(wektory_wlasne))
        tmp = []
        for i in czestosci_wlasne:
            tmp.append(wektory_wlasne[i[0]][1])
        return savetxt(str(self.lista_wektorow_q[0]) + '.', array(tmp).view(float))

"""
def start():
    return ZagadnienieWlasne(121, 1, 'FFT', 'II').wektory_wlasne()
"""


