from src.eig_problem.cProfiler import do_cprofile
from scipy.linalg import eig
from src.eig_problem.MacierzDoZagadnienia import MacierzDoZagadnienia
from src.eig_problem.ZagadnienieWlasne import ZagadnienieWlasne
import pandas as pd
import numpy as np
from math import *

class ZagadnienieWlasneAntidot(ZagadnienieWlasne):
    def __init__(self, ilosc_wektorow_q, skad_wspolczynnik='FFT'):

        ZagadnienieWlasne.__init__(self, ilosc_wektorow_q, skad_wspolczynnik='FFT')
        self.skad_wspolczynnik = skad_wspolczynnik

        self.tmp_table = pd.read_csv(self.input_fft, delimiter=' ', dtype=float, header=None).values
        re = self.tmp_table[:, 0::2]
        im = self.tmp_table[:, 1::2] * 1j
        self.wspolczynniki_fft = re + im

    def wspolczynniki_do_macierzy_material(self):
        # noinspection PyTypeChecker
        index1 = int(len(self.wspolczynniki_fft) / 2.) - int((self.ilosc_wektorow - 1) / 2)
        wspolczynniki = np.zeros((self.ilosc_wektorow, self.ilosc_wektorow), dtype=complex)
        for i in range(self.ilosc_wektorow):
            for j in range(self.ilosc_wektorow):
                wspolczynniki[i][j] = self.wspolczynniki_fft[i + index1][j + index1]
        assert len(wspolczynniki) == self.ilosc_wektorow
        assert wspolczynniki[int(self.ilosc_wektorow-1)/2][int(self.ilosc_wektorow - 1) / 2].imag == 0.
        return wspolczynniki

    def macierz_wspolczynnikow(self):
        macierz_antyprzekatna = self.wspolczynniki_do_macierzy_material
        # noinspection PyTypeChecker
        macierz_przekatna = np.zeros((len(macierz_antyprzekatna), len(macierz_antyprzekatna)))

        return np.concatenate((np.concatenate((macierz_antyprzekatna, macierz_przekatna), axis=0),
                    np.concatenate((macierz_przekatna, macierz_antyprzekatna), axis=0)), axis=1)

    def zagadnienie_wlasne(self, wektor_q, param):
        macierz_m = MacierzDoZagadnienia(self.skad_wspolczynnik).fill_matrix(wektor_q)
        #macierz_materialowa = self.macierz_wspolczynnikow()
        macierz_materialowa = np.identity(len(macierz_m)) # macierz identyczności - sprawdzenie poprawności metody
        return eig(macierz_m, macierz_materialowa, right=False, left=param)

    def czestosci_wlasne(self, wektor_q):
        assert type(wektor_q) == tuple, \
            'form of wektor_q is forbidden. wektor_q should be touple'
        assert len(wektor_q) == 2, \
            'form of wektor_q is forbidden. wektor_q should have two arguments'
        wartosci_wlasne = self.zagadnienie_wlasne(wektor_q, param=False)
        #wartosci_wlasne **= -1
        czestosci_wlasne = [i.imag * self.gamma * self.mu0H0 / 2.0 / pi for i in wartosci_wlasne if i.imag > 0]
        return list(sorted(czestosci_wlasne))[0:300]

    def wypisz_czestosci_do_pliku(self):
        plik = []
        for k in self.lista_wektorow_q:
            tmp = [k[0]]
            tmp.extend(self.czestosci_wlasne(k))
            plik.append(tmp)
        np.savetxt(str(self.a) + 'nm_' + str(self.ilosc_wektorow) + '_' + str(self.gamma/1e9) +'GHzid.txt', plik)

    def wektory_wlasne(self):
        """
        Metoda, której zadaniem jest wygenrowanie wektorów własnych, służących do wykreślenia profili wzbudzeń.
        :return: Plik txt zawierający639e38ab877de8ba05a68937e881e5316b0dcf86
        """
        assert len(self.lista_wektorow_q) == 1, 'Eigenvector should be calculated for only one position vector'
        wartosci_wlasne, wektory_wlasne = self.zagadnienie_wlasne(self.lista_wektorow_q[0], param=True)
        # noinspection PyTypeChecker
        wartosci_wlasne_index = np.argsort(wartosci_wlasne.imag**-1)
        wektory_wlasne = np.transpose(wektory_wlasne)
        wektory_wlasne = wektory_wlasne[wartosci_wlasne_index[self.ilosc_wektorow:]]
        return np.savetxt(str(self.lista_wektorow_q[0]) + '.', wektory_wlasne.view(float))

q = ZagadnienieWlasneAntidot(1)
#q.wektory_wlasne()
q.wypisz_czestosci_do_pliku()