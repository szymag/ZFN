import os
import glob
import re
import matplotlib.pyplot as plt
import numpy as np
from src.io.NumpyDataWriter import NumpyDataWriter
from src.fft_from_image.ElementaryUnitCellImage import ElementaryUnitCellImage


class FFTfromImage:
    def __init__(self):
        pass

    @staticmethod
    def fft(tablica):

        return np.fft.fftshift(np.fft.fft2(tablica, norm='ortho')) / len(tablica)

    def wywolaj_fft(self, path='.'):

        lista_plikow = ElementaryUnitCellImage(path).apply_for_every_image_in_directory()

        lista_fft = [self.fft(k) for k in lista_plikow]
        return lista_fft

    def wykres(self):
        lista_fft = self.wywolaj_fft
        for tablica in lista_fft():
            rozmiar = (len(tablica[0]), len(tablica))
            x, y = np.mgrid[slice(0, rozmiar[1], 1), slice(0, rozmiar[0], 1)]
            z = abs(tablica)
            plt.pcolor(x, y, z, cmap='gray')
            plt.colorbar()
            plt.show()

    def wypisz_do_pliku(self,data_writer,  path='', lista_fft=None):
        if lista_fft is None:
            lista_fft = self.wywolaj_fft(path)
        indeks = 1
        files = []
        for tablica in lista_fft:
            filepath = os.path.join(os.path.abspath(path),
                                    re.split(r'\.(?!\d)', str(list((glob.glob("*.png")))[indeks-1]))[0] + '.txt')
            files.append(filepath)
            np.savetxt(filepath, tablica.view(float))
            indeks += 1
        return files



FFTfromImage().wypisz_do_pliku()
