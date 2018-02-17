# -*- coding: utf-8 -*-
import os.path

import numpy as np
from scipy.linalg import eig
import matplotlib.pyplot as plt
from src.eig_problem.MacierzDoZagadnienia import MacierzDoZagadnienia

from src.eig_problem.ParametryMaterialowe import ParametryMaterialowe

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
        # self.mu0H0 = eval(mu0H0)
        self.mu0H0 = mu0H0
        self.H0 = self.mu0H0 / ParametryMaterialowe.mu0
        self.a = a
        self.lista_wektorow_q = [2 * np.pi * k / a for k in np.linspace(0.01, 0.99, ilosc_wektorow_q)]
        self.input_fft = os.path.join(scriptpath, input_fft)
        self.output_file = output_file
        # self.angle = eval(angle)
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
        # np.savetxt(self.output_file, self.czestosci_wlasne(self.lista_wektorow_q))
        plik = []
        for k in self.lista_wektorow_q:
            tmp = [k]
            tmp.extend(self.czestosci_wlasne(k))
            plik.append(tmp)
        # np.savetxt(self.output_file, plik)
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
    freq_vs_anlge = np.zeros((50, 91))
    for i in enumerate(np.arange(0, 91, 1)):
        ZagadnienieWlasne(1, 'heat_fft_Ni.txt', str(i[1]) + '.dat', angle=i[1]).wektory_wlasne()

        freq_vs_anlge[0:50, i[0]] = ZagadnienieWlasne(1, 'heat_fft_Ni.txt',
                                                      str(i[1]) + '.dat',
                                                      angle=i[1]).wypisz_czestosci_do_pliku()[0][1:]
        np.savetxt('freq_vs_angle_40.dat', freq_vs_anlge)


def start_1():
    fields = [0.030, 0.035, 0.040, 0.045, 0.05, 0.055, 0.06, 0.065]
    freq_vs_angle = np.zeros((50, len(fields)))

    m_min = 909
    for ind, field in enumerate(fields):
        ZagadnienieWlasne(1, 'c_coef_100.txt', str(field)
                          + '_' + str(m_min) + '.dat', mu0H0=field, angle=60).wektory_wlasne()
        freq_vs_angle[0:50, ind] = ZagadnienieWlasne(1, 'c_coef_100.txt',
                                                     'freq_vs_angle_' + str(field)
                                                     + '_' + str(m_min) + '.dat',
                                                     mu0H0=field,
                                                     angle=60).wypisz_czestosci_do_pliku()[0][1:]
    np.savetxt('freq_vs_angle_' + str(m_min) + '.dat', freq_vs_angle)


def start_2():
    fields = np.arange(0.015, 0.08, 0.005)
    angle = 10
    freq_vs_field = np.zeros((50, len(fields)))
    for field in enumerate(fields):
        freq_vs_field[0:50, field[0]] = ZagadnienieWlasne(1, 'c_coef_100.txt',
                                                          'freq_vs_field_' + str(angle) + '.dat',
                                                          mu0H0=field[1],
                                                          angle=angle).wypisz_czestosci_do_pliku()[0][1:]
    #np.savetxt('freq_vs_field_' + str(angle) + '.dat', freq_vs_field)
    for branch in range(50):
        plt.plot(fields, freq_vs_field[branch] / 1e9, color='blue')
    plt.xlabel('Field (T)')
    plt.ylabel('Freq (GHz)')
    plt.title('Modulation: M_max=.948, M_min=.793 (.5ns); sin-like profile')
    plt.show()

if __name__ == "__main__":
    start()
