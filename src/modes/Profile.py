# from src.eig_problem.cProfiler import do_cprofile
import ast
import glob
import os
import numpy as np
from src.eig_problem.ParametryMaterialowe import ParametryMaterialowe
from src.eig_problem.WektorySieciOdwrotnej import WektorySieciOdwrotnej
from cmath import *
import matplotlib.pyplot as plt


class Profile(ParametryMaterialowe):
    def __init__(self, ilosc_wektorow=49, typ_pola_wymiany=None, start_path="."):
        ParametryMaterialowe.__init__(self, ilosc_wektorow, typ_pola_wymiany)
        self.ilosc_wektorow = ilosc_wektorow
        self.lista_wektorow = WektorySieciOdwrotnej(self.a, self.b, ilosc_wektorow).lista_wektorow('min')
        self.sciezka = glob.glob(os.path.join(start_path, "*."))
        self.wektory_wlasne = np.loadtxt(self.sciezka[0]).view(complex)
        self.wektor_q = ast.literal_eval(self.sciezka[0].strip('.')[1:])

    def magnetyzacja_w_punkcie(self, wektor_r, numer_modu):
        assert len(wektor_r) == 2
        lista_wektorow = enumerate(self.lista_wektorow)
        tmp = 0.
        # TODO: Sprawdzić porządek wektorów sieci odwr poprzez odtworzenie parametrów materiałowych
        for wektor_odwr in lista_wektorow:
            tmp += self.wektory_wlasne[numer_modu-1][wektor_odwr[0]] * exp(1.j * ((wektor_odwr[1][0] +
                    self.wektor_q[0]) * wektor_r[0] + (wektor_odwr[1][1] + self.wektor_q[1]) * wektor_r[1]))
        return abs(tmp)

    def mapa_profile(self, numer_modu):
        print(self.lista_wektorow)
        lista_x = np.linspace(0, self.a*2, 200)
        lista_y = np.linspace(0, self.b*2, 200)
        lista_wartosci = np.zeros((len(lista_x), len(lista_y)))
        for x in enumerate(lista_x):
            for y in enumerate(lista_y):
                lista_wartosci[x[0]][y[0]] = self.magnetyzacja_w_punkcie((x[1], y[1]), numer_modu)
        return lista_x, lista_y, lista_wartosci

    def wykreslenie_profili(self, numer_modu):
        lista_x, lista_y, lista_wartosci = self.mapa_profile(numer_modu)
        x, y = np.meshgrid(np.array(lista_x), np.array(lista_y))
        plt.pcolor(x, y, np.array(lista_wartosci))
        plt.colorbar()
        plt.show()

print(Profile().wykreslenie_profili(1))
