import os
import glob
import re
import matplotlib.pyplot as plt
import numpy as np

from src.fft_from_image.TablicaWartosciPikseli import TablicaWartosciPikseli


class FFT:
    """
    Klasa, której zadaniem jest przygotowanie plików '*.txt' zawierających wartości współczynników Fouriera dla zadanych
    plików graficznych.
    """

    def __init__(self):
        pass

    @staticmethod
    def fft(tablica):
        """
        Wywołanie FFT jako metody obliczającej współczynniki Fouriera dla zadanej tablicy.
        :param tablica: Tablica wartości. Standardowo, jest to tablica wartości pikseli.
        :return: Wspołczynniki Fouriera w postaci tablicy, o wymiarach takich jak zadana tablica.
        """

        return np.fft.fftshift(np.fft.fft2(tablica, norm='ortho')) / len(tablica)

    def wywolaj_fft(self, path='.'):
        """
        Dla każdego pliku '*.png' znajdującego się w tym samym katalogu co klasa, wywoływana jest metdoda fft
        :param path: Scieżka pliku
        :return: Lista tablic, które pochodzą z wywołania metody fft.
        """
        lista_plikow = TablicaWartosciPikseli(path).tablica_dla_plikow()

        lista_fft = [self.fft(k) for k in lista_plikow]
        return lista_fft

    def wykres(self):
        """
        W celu sprawdzenia poprawności, wykreślane są wykresie typu 'density plot' amplitudy współczynników Fouriera.
        :return: Dla każdego pliku w katalogu '.*png' wyświetlony będzie wykres.
        """
        lista_fft = self.wywolaj_fft
        for tablica in lista_fft():
            rozmiar = (len(tablica[0]), len(tablica))
            x, y = np.mgrid[slice(0, rozmiar[1], 1), slice(0, rozmiar[0], 1)]
            z = abs(tablica)
            plt.pcolor(x, y, z, cmap='gray')
            plt.colorbar()
            plt.show()

    def wypisz_do_pliku(self, path='', lista_fft=None):
        """
        Metoda wypisująca do plików w postaci tablic, współczynniki Fouriera.
        Odpowiednio dla każdego pliku '*.png' w katalogu.
        :param lista_fft:
        :param path: Ścieżka
        :return: Pliki typu '*.txt' dla każdego obrazka '*.png'
        """
        if lista_fft is None:
            lista_fft = self.wywolaj_fft(path)
        indeks = 1
        files = []
        for tablica in lista_fft:
            # TODO: usprawnić nazywanie plików
            filepath = os.path.join(os.path.abspath(path),
                                    re.split(r'\.(?!\d)', str(list((glob.glob("*.png")))[indeks-1]))[0] + '.txt')
            files.append(filepath)
            np.savetxt(filepath, tablica.view(float))
            indeks += 1
        return files


FFT().wypisz_do_pliku()
