# -*- coding: utf-8 -*-
import os.path

import numpy as np
from scipy.linalg import eig
from src.eig_problem.MacierzDoZagadnienia import MacierzDoZagadnienia

from src.modes.ParametryMaterialowe import ParametryMaterialowe

scriptpath = os.path.dirname(__file__)

class ZagadnienieWlasne:


    def __init__(self, ilosc_wektorow_q, input_fft, output_file, mu0H0=ParametryMaterialowe.mu0H0,
                  a=ParametryMaterialowe.a, gamma=ParametryMaterialowe.gamma, angle=ParametryMaterialowe.angle):
        """
        :param ilosc_wektorow_q: Odpowiada za gęstość siatki, na wykresie dyspersji.
        :skad_wspolczynnik: Argument, który odpowiada z źródło pochodzenia wspoółczynników Fouriera. Możliwe wartości
        to DFT oraz FFT.
        """
        self.gamma = gamma
        #self.mu0H0 = eval(mu0H0)
        self.mu0H0 = mu0H0
        self.H0 = self.mu0H0 / ParametryMaterialowe.mu0
        self.a = a
        self.lista_wektorow_q = [2 * np.pi * k / a for k in np.linspace(0.0001, 0.99, ilosc_wektorow_q)]
        self.input_fft = os.path.join(scriptpath, input_fft)
        self.output_file = output_file
        #self.angle = eval(angle)
        self.angle = angle

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
        macierz_m = MacierzDoZagadnienia(self.input_fft, wektor_q,
                                         angle=self.angle, H0=self.H0).matrix_angle_dependence(wektor_q)

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
        return plik

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

    #return ZagadnienieWlasne(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], angle=sys.argv[5]).wypisz_czestosci_do_pliku()
    ZagadnienieWlasne(30, 'c_coef_100.txt', 'dys_090.dat').wypisz_czestosci_do_pliku()
    #for i in range(0, 92, 2):
    #freq_vs_anlge = np.zeros((50, 90))
    #for i in range(0, 90, 1):
    #    ZagadnienieWlasne(1, 'heat_fft_t700.txt', 't700_' + str(i) + '.dat', angle=i).wektory_wlasne()
    #    #print(len(ZagadnienieWlasne(1, 'heat_fft.txt', 'fmr'+str(i)+'.dat', angle=i).wypisz_czestosci_do_pliku()[0][1:]))
    #    freq_vs_anlge[0:50, i] = ZagadnienieWlasne(1, 'heat_fft_t700.txt', 'fmr'+str(i)+'.dat', angle=i).wypisz_czestosci_do_pliku()[0][1:]
    #    np.savetxt('freq_vs_angle_t700.dat', freq_vs_anlge)

if __name__ == "__main__":
    start()
