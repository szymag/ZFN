# from src.eig_problem.cProfiler import do_cprofile
from scipy.linalg import eig
from src.eig_problem.MacierzDoZagadnienia import MacierzDoZagadnienia
from src.eig_problem.ZagadnienieWlasne import ZagadnienieWlasne
import numpy as np


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
        """
        Metoda tworząca listę współczynników. Na podstawie położenia w tablicy, określane jest położenie w liście
        :return: Lista współczynników.
        """
        index1 = int(len(self.wspolczynniki_fft) / 2.) - self.ilosc_wektorow
        wspolczynniki = np.zeros((self.ilosc_wektorow, self.ilosc_wektorow), dtype=complex)
        for i in range(self.ilosc_wektorow):
            for j in range(self.ilosc_wektorow):
                wspolczynniki[i][j] = self.wspolczynniki_fft[i + index1][j + index1]
        assert len(wspolczynniki) == self.ilosc_wektorow
        return wspolczynniki

    def macierz_material(self):
        macierz_antyprzekatna = self.wspolczynniki_do_macierzy_material()
        macierz_przekatna = np.zeros((len(macierz_antyprzekatna), len(macierz_antyprzekatna)))
        return np.concatenate((np.concatenate((macierz_przekatna, macierz_antyprzekatna), axis=1),
                               np.concatenate((macierz_antyprzekatna, macierz_przekatna), axis=1)), axis=0)

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
        macierz_materialowa = self.macierz_material()
        return eig(macierz_materialowa, macierz_m, right=param)


