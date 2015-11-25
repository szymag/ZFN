import glob
import os

from PIL import Image
from numpy import zeros


class TablicaWartosciPikseli:
    def __init__(self, start_path="."):
        self.lista_plikow = list((glob.glob(os.path.join(start_path, "*.png"))))

    @staticmethod
    def stworz_tablice(wczytany_plik):
        plik = Image.open(wczytany_plik)
        piksele = [abs((k / 255) - 1) + 1 for k in list(plik.getdata(0))]
        rozmiar = plik.size
        tablica_pikseli = list(zeros(rozmiar[1]))
        for i in range(0, rozmiar[1]):
            tablica_pikseli[i] = piksele[rozmiar[0] * i:(i + 1) * rozmiar[0]]
        return tablica_pikseli

    def tablica_dla_plikow(self):
        return [self.stworz_tablice(str(k)) for k in self.lista_plikow]

    def konwersja_dwa_kolory(self):
        return self.lista_plikow
