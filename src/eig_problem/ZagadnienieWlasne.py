# -*- coding: utf-8 -*-
import numpy as np
from scipy.linalg import eig
import sys
from src.eig_problem.MacierzDoZagadnienia import MacierzDoZagadnienia
from src.eig_problem.ParametryMaterialowe import ParametryMaterialowe
from src.eig_problem.cProfiler import do_cprofile
import os.path

scriptpath = os.path.dirname(__file__)

class ZagadnienieWlasne:


    def __init__(self, ilosc_wektorow_q, input_fft, output_file, mu0H0=ParametryMaterialowe.mu0H0,
                  a=ParametryMaterialowe.a, gamma=ParametryMaterialowe.gamma):
        """
        :param ilosc_wektorow_q: Odpowiada za gęstość siatki, na wykresie dyspersji.
        :skad_wspolczynnik: Argument, który odpowiada z źródło pochodzenia wspoółczynników Fouriera. Możliwe wartości
        to DFT oraz FFT.
        """
        self.gamma = gamma
        self.mu0H0 = eval(mu0H0)
        #self.mu0H0 = mu0H0
        self.H0 = self.mu0H0 / ParametryMaterialowe.mu0
        self.a = a
        self.lista_wektorow_q = [2 * np.pi * k / a for k in np.linspace(0.0001, 0.999, ilosc_wektorow_q)]
        self.input_fft = os.path.join(scriptpath, input_fft)
        self.output_file = output_file

    # @do_cprofile
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
        macierz_m = MacierzDoZagadnienia(self.input_fft, wektor_q, self.H0).matrix_angle_dependence(wektor_q)
        #macierz_m = MacierzDoZagadnienia(self.input_fft, wektor_q, H0=self.H0).matrix_gen_backward_volume(wektor_q)

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
        wartosci_wlasne = self.zagadnienie_wlasne(wektor_q, param=False)
        czestosci_wlasne = [i.imag * self.gamma * self.mu0H0 / 2.0 / np.pi for i in wartosci_wlasne if i.imag > 0]

        return list(sorted(czestosci_wlasne)[:50])

    def wypisz_czestosci_do_pliku(self):
        """
        Metoda, która wypisuje w kolumnach wyniki. W pierwszej jest wektor blochowski, a w dalszych odpowiadające mu
        częstotliwości własne.
        :return: Plik tekstowy.
        """
        #np.savetxt(self.output_file, self.czestosci_wlasne(self.lista_wektorow_q))
        plik = []
        for k in self.lista_wektorow_q:
            tmp = [k]
            tmp.extend(self.czestosci_wlasne(k))
            plik.append(tmp)
        np.savetxt(self.output_file, plik)

    def wektory_wlasne(self):
        """
        Metoda, której zadaniem jest wygenrowanie wektorów własnych, służących do wykreślenia profili wzbudzeń.
        :return: Plik txt zawierający
        """
        # assert len(self.lista_wektorow_q) == 1, 'Eigenvector should be calculated for only one position vector'
        wartosci_wlasne, wektory_wlasne = self.zagadnienie_wlasne(self.lista_wektorow_q[0], True)
        wartosci_wlasne_index = np.argsort(wartosci_wlasne.imag)
        wektory_wlasne = np.transpose(wektory_wlasne)
        wektory_wlasne = wektory_wlasne[wartosci_wlasne_index[sum(wartosci_wlasne.imag < 0):]]
        np.savetxt(self.output_file, wektory_wlasne.view(float))
        return wektory_wlasne


def start():

     return ZagadnienieWlasne(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]).wypisz_czestosci_do_pliku()
    # return ZagadnienieWlasne(30, 'p_coef_100*2.txt', 'dys_90.dat').wypisz_czestosci_do_pliku()
    # return ZagadnienieWlasne(1, 'p_coef_100*2.txt', 'dys_0.dat').wektory_wlasne()

if __name__ == "__main__":
    start()
