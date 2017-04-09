# -*- coding: utf-8 -*-
import numpy as np
from scipy.linalg import eig
import sys
from src.eig_problem.EigenMatrix1D import EigenMatrix1D
from src.eig_problem.InputParameter import InputParameter
import os.path

scriptpath = os.path.dirname(__file__)

class ZagadnienieWlasne:


    def __init__(self, ilosc_wektorow_q, input_fft, output_file, mu0H0=InputParameter.mu0H0,
                  a=InputParameter.a, gamma=InputParameter.gamma, angle=InputParameter.angle):
        self.gamma = gamma
        #self.mu0H0 = eval(mu0H0)
        self.mu0H0 = mu0H0
        #self.H0 = self.mu0H0 / ParametryMaterialowe.mu0
        self.a = a
        self.lista_wektorow_q = [2 * np.pi * k / a for k in np.linspace(0.0001, 0.9999, ilosc_wektorow_q)]
        self.input_fft = os.path.join(scriptpath, input_fft)
        self.output_file = output_file
        self.angle = angle

    # @do_cprofile
    def zagadnienie_wlasne(self, wektor_q, param):
        macierz_m = EigenMatrix1D(self.input_fft, wektor_q,
                                         angle=self.angle).matrix_angle_dependence(wektor_q)

        return eig(macierz_m, right=param)  # trzeba pamiętać o włączeniu/wyłączeniu generowania wektorów

    def czestosci_wlasne(self, wektor_q):
        wartosci_wlasne = self.zagadnienie_wlasne(wektor_q, param=False)
        czestosci_wlasne = [i.imag * self.gamma * self.mu0H0 / 2.0 / np.pi for i in wartosci_wlasne if i.imag > 0]

        return list(sorted(czestosci_wlasne)[:50])

    def wypisz_czestosci_do_pliku(self):
        #np.savetxt(self.output_file, self.czestosci_wlasne(self.lista_wektorow_q))
        plik = []
        for k in self.lista_wektorow_q:
            tmp = [k]
            tmp.extend(self.czestosci_wlasne(k))
            plik.append(tmp)
        np.savetxt(self.output_file, plik)

    def wektory_wlasne(self):
        # assert len(self.lista_wektorow_q) == 1, 'Eigenvector should be calculated for only one position vector'
        wartosci_wlasne, wektory_wlasne = self.zagadnienie_wlasne(self.lista_wektorow_q[0], True)
        wartosci_wlasne_index = np.argsort(wartosci_wlasne.imag)
        wektory_wlasne = np.transpose(wektory_wlasne)
        wektory_wlasne = wektory_wlasne[wartosci_wlasne_index[sum(wartosci_wlasne.imag < 0):]]
        np.savetxt(self.output_file, wektory_wlasne.view(float))
        return wektory_wlasne


def start():

    # return ZagadnienieWlasne(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]).wypisz_czestosci_do_pliku()
    # return ZagadnienieWlasne(30, 'c_coef_100.txt', 'dys_90.dat').wypisz_czestosci_do_pliku()
    for i in range(0, 93, 3):
        ZagadnienieWlasne(1, 'c_coef_100.txt', 'vec_' + str(i) + '.dat', angle=i).wektory_wlasne()


if __name__ == "__main__":
    start()
