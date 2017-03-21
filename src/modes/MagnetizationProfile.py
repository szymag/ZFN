# from src.eig_problem.cProfiler import do_cprofile
import ast
import glob
import os

import matplotlib.pyplot as plt
import numpy as np

from src.eig_problem.ParametryMaterialowe import ParametryMaterialowe
from src.eig_problem.WektorySieciOdwrotnej import WektorySieciOdwrotnej


class Profile2D:
    def __init__(self, ilosc_wektorow=441, start_path="."):
        self.ilosc_wektorow = ilosc_wektorow
        self.lista_wektorow = WektorySieciOdwrotnej(self.a, self.b, ilosc_wektorow).lista_wektorow('min')
        self.sciezka = glob.glob(os.path.join(start_path, "*."))
        self.wektory_wlasne = np.loadtxt(self.sciezka[0]).view(complex)
        self.wektor_q = ast.literal_eval(self.sciezka[0].strip('.')[1:])

    def magnetyzacja_w_punkcie(self, wektor_r, numer_modu):
        wektory_wlasne = np.array(self.wektory_wlasne[numer_modu-1])
        wektory_odwrotne = np.array(self.lista_wektorow)
        return abs(np.sum(wektory_wlasne[0:self.ilosc_wektorow] *
                          np.prod(np.exp(1j * wektory_odwrotne * (wektor_r + self.wektor_q)), axis=1)))

    def mapa_profile(self, numer_modu, dokladnosc):
        x = np.linspace(0, self.a, dokladnosc)
        y = np.linspace(0, self.a, dokladnosc)
        m = np.zeros(dokladnosc * dokladnosc)
        for i in enumerate(np.dstack(np.meshgrid(x, y)).reshape(-1, 2)):
            m[i[0]] = self.magnetyzacja_w_punkcie(i[1], numer_modu)
        return x, y, m.reshape((dokladnosc, dokladnosc))

    def wykreslenie_profili(self, numer_modu, dokladnosc):
        lista_x, lista_y, lista_wartosci = self.mapa_profile(numer_modu, dokladnosc)
        x, y = np.meshgrid(np.array(lista_x), np.array(lista_y))
        plt.pcolor(x, y, np.array(lista_wartosci))
        plt.colorbar()
        plt.show()

# Profile2D().wykreslenie_profili(1, 150)

class Profile1D:
    def __init__(self, mode_number):
        self.eig_vectors = np.loadtxt('dys_90.dat').view(complex)
        self.lattice_const = ParametryMaterialowe.a
        self.mode_number = mode_number - 1

    def plot(self):
        tmp = self.spatial_distribution_dynamic_magnetization(500)
        #tmp1 = self.elementary_cell_reconstruction(500)
        plt.plot(tmp[0], tmp[1])
        #plt.plot(tmp1[0], tmp1[1])
        #plt.ylim((0,3.5))
        plt.savefig('90deg_mode_' + str(self.mode_number) + '.svg')
        plt.clf()

    def spatial_distribution_dynamic_magnetization(self, grid):
        mode = self.eig_vectors[self.mode_number, 0:100]
        x = np.linspace(-self.lattice_const, 0, grid)
        tmp = np.zeros(grid)
        for i in enumerate(x):
            tmp[i[0]] = self.inverse_discrete_foutier_transform(mode, i[1])
        return x, tmp

    def elementary_cell_reconstruction(self, grid):
        coefficient = np.transpose(np.loadtxt('p_coef_100*2.txt').view(complex))
        x = np.linspace(-self.lattice_const, 0, grid)
        tmp = np.zeros(grid)
        for i in enumerate(x):
            tmp[i[0]] = self.inverse_discrete_foutier_transform(coefficient, i[1])
        return x, (tmp+.3) /1.3

    def inverse_discrete_foutier_transform(self, data, vector_position):
        reciprocal_vectors = np.array(2 * np.pi * WektorySieciOdwrotnej(max(data.shape)).lista_wektorow1d('min')
                                      / self.lattice_const)
        return abs(np.sum(data * np.exp(1j * reciprocal_vectors * vector_position))) **2

for i in range(1, 5):
    Profile1D(i).plot()