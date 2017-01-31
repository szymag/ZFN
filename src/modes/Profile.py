# from src.eig_problem.cProfiler import do_cprofile
import ast, glob, os
import numpy as np
from src.eig_problem.ParametryMaterialowe import ParametryMaterialowe
from src.eig_problem.WektorySieciOdwrotnej import WektorySieciOdwrotnej
import matplotlib.pyplot as plt


class Profile(ParametryMaterialowe):
    def __init__(self, ilosc_wektorow=441, start_path="."):
        ParametryMaterialowe.__init__(self)
        self.ilosc_wektorow = ilosc_wektorow
        self.lista_wektorow = WektorySieciOdwrotnej(self.a, self.b, ilosc_wektorow).lista_wektorow2d('min')
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

Profile().wykreslenie_profili(1, 150)
