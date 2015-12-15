import ast
import glob
import os
import numpy as np
from src.eig_problem.ParametryMaterialowe import ParametryMaterialowe
from src.eig_problem.WektorySieciOdwrotnej import WektorySieciOdwrotnej
from math import *

class Profile(ParametryMaterialowe):
    def __init__(self, ilosc_wektorow=25, typ_pola_wymiany=None, start_path="."):
        ParametryMaterialowe.__init__(self, ilosc_wektorow, typ_pola_wymiany)
        self.ilosc_wektorow = ilosc_wektorow
        self.lista_wektorow = WektorySieciOdwrotnej(self.a, self.b, ilosc_wektorow).lista_wektorow('min')
        self.sciezka = glob.glob(os.path.join(start_path, "*."))
        self.plik = np.loadtxt(self.sciezka[0]).view(complex)
        self.wektor_q = ast.literal_eval(self.sciezka[0].strip('.')[1:])

    def wspolrzedna(self, wektor_r, numer_modu):
        assert len(wektor_r) == 2
        lista_wektorow = enumerate(self.lista_wektorow)
        tmp = 0.
        for i in lista_wektorow:
            pass
        pass

    def q(self):
        print(len(self.plik[0]))


Profile().q()
