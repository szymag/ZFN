import os
import matplotlib.pyplot as plt
import numpy

from src.fft_from_image.TablicaWartosciPikseli import TablicaWartosciPikseli


class FFT:
    def __init__(self):
        pass

    def fft(self, tablica):
        return numpy.fft.fftshift(numpy.fft.fft2(tablica))

    def wywolaj_fft(self, path='.'):
        lista_plikow = TablicaWartosciPikseli(path).tablica_dla_plikow()
        lista_fft = [self.fft(k) for k in lista_plikow]
        return lista_fft

    def wykres(self):
        lista_fft = self.wywolaj_fft
        for tablica in lista_fft():
            rozmiar = (len(tablica[0]), len(tablica))
            x, y = numpy.mgrid[slice(0, rozmiar[1], 1), slice(0, rozmiar[0], 1)]
            z = tablica
            plt.pcolor(x, y, z, cmap='gray')
            plt.colorbar()
            plt.show()

    def wypisz_do_pliku(self, path='', lista_fft=None):
        if lista_fft is None:
            lista_fft = self.wywolaj_fft(path)
        indeks = 1
        files = []
        for tablica in lista_fft:
            # TODO: usprawnić nazywanie plików
            filepath = os.path.join(os.path.abspath(path), 'fft' + str(indeks) +'.txt')
            files.append(filepath)
            numpy.savetxt(filepath, tablica.view(float))
            indeks += 1
        return files

    def wczytaj_z_pliku(self):
        tablica = numpy.loadtxt(".txt")
        return tablica
