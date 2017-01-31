from src.eig_problem.cProfiler import do_cprofile
import numpy as np
from math import sqrt
from scipy.linalg import eig
from src.eig_problem.MacierzDoZagadnienia import MacierzDoZagadnienia
from src.eig_problem.ParametryMaterialowe import ParametryMaterialowe
import sys


class ZagadnienieWlasne:
    """
    Klasa, w której zdefiniowane jest uogólnione zagadnienie własne. Dziedziczy ona po 'Parametry materiałowe, gdyż
    potrzebne są w niej infoemacje o strukturze kryształu magnonicznego.
    """

    def __init__(self, ilosc_wektorow_q, direction, a=ParametryMaterialowe.a, b=ParametryMaterialowe.b,
                 gamma=ParametryMaterialowe.gamma, mu0H0=ParametryMaterialowe.mu0H0):
        """
        :param ilosc_wektorow_q: Odpowiada za gęstość siatki, na wykresie dyspersji.
        """
        self.number_of_vec_q = ilosc_wektorow_q
        self.gamma = gamma
        self.mu0H0 = mu0H0
        self.direction = direction
        self.a = a
        self.b = b
        self.start_vec_q = 0.01
        self.end_vec_q = 0.002

    def list_vector_q(self):
        if self.direction == 'y':
            return [[(2 * np.pi * k / self.a), 0]
                    for k in np.linspace(self.start_vec_q, self.end_vec_q, self.number_of_vec_q)]
        elif self.direction == 'x':
            return [[0, (2 * np.pi * k / self.b)]
                    for k in np.linspace(self.start_vec_q, self.end_vec_q, self.number_of_vec_q)]
        elif self.direction == 'xy':
            return [[(2 * np.pi * k / self.a), (2 * np.pi * k / self.b)]
                    for k in np.linspace(self.start_vec_q, self.end_vec_q, self.number_of_vec_q)]
        else:
            sys.exit('Wrong argumnet for direction was set')

#    @do_cprofile
    def zagadnienie_wlasne(self, wektor_q, param):
        """
        Metoda, która wywołuje algorytm rozwiązywania zagadnienia własnego. Tworzy sobie tablicę,
        dla której następnie oblicza wartości i wektory własne. Jest ona także przystosowana dla
        uogólnionego zagadnienia własnego.
        :param param:
        :type wektor_q tuple
        :param wektor_q: Blochowski wektor. Jest on "uciąglony". Jest on zmienną przy wyznaczaniu dyspersji.
        :return: Wartości własne. Wektory własne są obecnie wyłączone.
        """
        macierz_m = MacierzDoZagadnienia("penrose.txt", 441, wektor_q).wypelnienie_macierzy()
        return eig(macierz_m, right=param)  # trzeba pamiętać o włączeniu/wyłączeniu generowania wektorów

    def czestosci_wlasne(self, wektor_q):
        """
        Metoda, w której przelicza się wartości własne w metodzie 'zagadnienie wlasne' na czestosci wlasne' wzbudzen
        w krysztale magnonicznym.
        :param wektor_q:  Blochowski wektor. Jest on "uciąglony". Jest on zmienną przy wyznaczaniu dyspersji.
        :return: Częstości własne, obliczane według wzoru: w = -a.imag * gamma * mu0*H0, gdzie w - cżęstość własne,
        a.imag część urojona wartości własnej. Ponieważ wartości własne są przmnażane przez 1/i, wystarczy wziąć część
        zespoloną, by uzysać częstotliwość rzeczywistą. Część zespolona w częstotliwośći odpowiada za tłumienie
        w układzie, a takiego się nie zakłada, więc powinny one być stosunkowo niewielkie - pomijalne.
        """
        assert len(wektor_q) == 2, \
            'form of wektor_q is forbidden. wektor_q should have two arguments'
        wartosci_wlasne = self.zagadnienie_wlasne(wektor_q, param=False)
        czestosci_wlasne = [i.imag * self.gamma * self.mu0H0 / 2.0 / np.pi for i in wartosci_wlasne if i.imag > 0]

        return list(sorted(czestosci_wlasne)[:50])

    def wypisz_czestosci_do_pliku(self):
        """
        Metoda, która wypisuje w kolumnach wyniki. W pierwszej jest wektor blochowski, a w dalszych odpowiadające mu
        częstotliwości własne.
        :return: Plik tekstowy.
        """
        plik = []
        if self.direction == 'x':
            for k in self.list_vector_q():
                tmp = [k[1]]
                tmp.extend(self.czestosci_wlasne(k))
                plik.append(tmp)
        elif self.direction == 'y':
            for k in self.list_vector_q():
                tmp = [k[0]]
                tmp.extend(self.czestosci_wlasne(k))
                plik.append(tmp)
        elif self.direction == 'xy':
            for k in self.list_vector_q():
                tmp = [sqrt(k[0] ** 2 + k[1] ** 2)]
                tmp.extend(self.czestosci_wlasne(k))
                plik.append(tmp)
        np.savetxt('test.txt', plik)

    def wektory_wlasne(self):
        """
        Metoda, której zadaniem jest wygenrowanie wektorów własnych, służących do wykreślenia profili wzbudzeń.
        :return: Plik txt zawierający
        """
        assert len(self.list_vector_q()) == 1, 'Eigenvector should be calculated for only one position vector'
        wartosci_wlasne, wektory_wlasne = self.zagadnienie_wlasne(self.list_vector_q()[0], param=True)
        wartosci_wlasne_index = np.argsort(wartosci_wlasne.imag)
        wektory_wlasne = np.transpose(wektory_wlasne)
        wektory_wlasne = wektory_wlasne[wartosci_wlasne_index[self.ilosc_wektorow:]]
        return np.savetxt(str(self.list_vector_q()[0]) + '.', wektory_wlasne.view(float))


def start():
    return ZagadnienieWlasne(1, 'y').wypisz_czestosci_do_pliku()

if __name__ == "__main__":
    start()
