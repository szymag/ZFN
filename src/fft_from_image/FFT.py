import matplotlib.pyplot as plt
from numpy import fft, mgrid, savetxt, loadtxt

from src.fft_from_image.TablicaWartosciPikseli import TablicaWartosciPikseli


class FFT(TablicaWartosciPikseli):
    def __init__(self):
        TablicaWartosciPikseli.__init__(self)

    @staticmethod
    def fft(tablica):
        return fft.fftshift(abs(fft.fft2(tablica, norm="ortho")))

    def wywolaj_fft(self):
        lista_plikow = TablicaWartosciPikseli().tablica_dla_plikow()
        lista_fft = [self.fft(k) for k in lista_plikow]
        return lista_fft

    def wykres(self):
        lista_fft = self.wywolaj_fft()
        for tablica in lista_fft:
            rozmiar = (len(tablica[0]), len(tablica))
            x, y = mgrid[slice(0, rozmiar[1], 1), slice(0, rozmiar[0], 1)]
            z = tablica
            plt.pcolor(x, y, z, cmap='gray')
            plt.colorbar()
            plt.show()

    def wypisz_do_pliku(self):
        lista_fft = self.wywolaj_fft()
        indeks = 1
        for tablica in lista_fft:
            # TODO: usprawnić nazywanie plików
            savetxt(str(indeks) + ".txt", tablica)
            indeks += 1

    def wczytaj_z_pliku(self):
        tablica = loadtxt("1.txt")
        return tablica
