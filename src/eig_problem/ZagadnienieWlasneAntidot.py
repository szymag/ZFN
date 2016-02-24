from src.eig_problem.cProfiler import do_cprofile
from scipy.linalg import eig
from src.eig_problem.MacierzDoZagadnienia import MacierzDoZagadnienia
from src.eig_problem.ZagadnienieWlasne import ZagadnienieWlasne
import numpy as np
from math import *


class ZagadnienieWlasneAntidot(ZagadnienieWlasne):
    def __init__(self, ilosc_wektorow, ilosc_wektorow_q, skad_wspolczynnik='FFT', typ_pole_wymiany='II',
                 filepath='fft1.txt'):

        ZagadnienieWlasne.__init__(self, ilosc_wektorow, ilosc_wektorow_q, skad_wspolczynnik='FFT',
                                   typ_pole_wymiany='II')
        self.typ_pole_wymiany = typ_pole_wymiany
        self.skad_wspolczynnik = skad_wspolczynnik
        self.wspolczynniki_fft = np.loadtxt(filepath).view(complex)
        self.ilosc_wektorow = ilosc_wektorow

    def wspolczynniki_do_macierzy_material(self):
        index1 = int(len(self.wspolczynniki_fft) / 2.) - int((self.ilosc_wektorow - 1)/2)
        wspolczynniki = np.zeros((self.ilosc_wektorow, self.ilosc_wektorow), dtype=complex)
        for i in range(self.ilosc_wektorow):
            for j in range(self.ilosc_wektorow):
                wspolczynniki[i][j] = self.wspolczynniki_fft[i + index1][j + index1]
        assert len(wspolczynniki) == self.ilosc_wektorow
        assert wspolczynniki[int(self.ilosc_wektorow-1)/2][int(self.ilosc_wektorow-1)/2].imag == 0.
        return wspolczynniki

    def macierz_material(self):
        macierz_antyprzekatna = self.wspolczynniki_do_macierzy_material()
        macierz_przekatna = np.zeros((len(macierz_antyprzekatna), len(macierz_antyprzekatna)))

        return np.concatenate((np.concatenate((macierz_antyprzekatna, macierz_przekatna), axis=0),
                               np.concatenate((macierz_przekatna, macierz_antyprzekatna), axis=0)), axis=1)

    def zagadnienie_wlasne(self, wektor_q, param):
        macierz_m = MacierzDoZagadnienia(self.ilosc_wektorow, self.skad_wspolczynnik,
                                         self.typ_pola_wymiany).wypelnienie_macierzy(wektor_q)
        #macierz_materialowa = self.macierz_material()
        macierz_materialowa = np.identity(len(macierz_m)) # macierz identyczności - sprawdzenie poprawności metody
        return eig(macierz_m, macierz_materialowa, right=param)

    def czestosci_wlasne(self, wektor_q):
        assert type(wektor_q) == tuple, \
            'form of wektor_q is forbidden. wektor_q should be touple'
        assert len(wektor_q) == 2, \
            'form of wektor_q is forbidden. wektor_q should have two arguments'
        wartosci_wlasne = self.zagadnienie_wlasne(wektor_q, param=False)
        # wartosci_wlasne, wektory_wlasny = self.zagadnienie_wlasne(wektor_q)
        # savetxt(str(wektor_q), wektory_wlasny.view(float))
        czestosci_wlasne = [i.imag * self.gamma * self.mu0H0 / 2.0 / pi for i in wartosci_wlasne if i.imag > 0]
        return list(sorted(czestosci_wlasne))[0:400]

    def wypisz_czestosci_do_pliku(self):
        plik = []
        for k in self.lista_wektorow_q:
            tmp = [k[0]]
            tmp.extend(self.czestosci_wlasne(k))
            plik.append(tmp)
        np.savetxt('A1_9xx_176.txt', plik)


q = ZagadnienieWlasneAntidot(841, 50)
#print(q.macierz_material())
q.wypisz_czestosci_do_pliku()