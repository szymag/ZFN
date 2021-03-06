import glob
import os
import re

import matplotlib.pyplot as plt
import numpy as np

from src.fft_from_image.Sequences import ThueMorse, Fibonacci, Periodic, Heated, Custom
from src.fft_from_image.ElementaryUnitCellImage import ElementaryUnitCellImage



class FFT:
    """
    Klasa, której zadaniem jest przygotowanie plików '*.txt' zawierających wartości współczynników Fouriera dla zadanych
    plików graficznych.
    """

    def __init__(self):
        pass

    @staticmethod
    def fft2d(tablica):
        """
        Wywołanie FFT jako metody obliczającej współczynniki Fouriera dla zadanej tablicy.
        :param tablica: Tablica wartości. Standardowo, jest to tablica wartości pikseli.
        :return: Wspołczynniki Fouriera w postaci tablicy, o wymiarach takich jak zadana tablica.
        """
        return np.fft.fftshift(np.fft.fft2(tablica, norm='ortho')) / len(tablica)

    def wywolaj_fft2d(self, path='.'):
        """
        Dla każdego pliku '*.png' znajdującego się w tym samym katalogu co klasa, wywoływana jest metdoda fft
        :param path: Scieżka pliku

        :return: Lista tablic, które pochodzą z wywołania metody fft.
        """
        lista_plikow = ElementaryUnitCellImage(path).apply_for_every_image_in_directory()
        lista_fft = [self.fft2d(k) for k in lista_plikow]
        return lista_fft

    def wykres(self):
        """
        W celu sprawdzenia poprawności, wykreślane są wykresie typu 'density plot' amplitudy współczynników Fouriera.
        :return: Dla każdego pliku w katalogu '.*png' wyświetlony będzie wykres.
        """
        lista_fft = self.wywolaj_fft2d
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
            lista_fft = self.wywolaj_fft2d(path)
        indeks = 1
        files = []
        for tablica in lista_fft:
            filepath = os.path.join(os.path.abspath(path),
                                    re.split(r'\.(?!\d)', str(list((glob.glob("*.png")))[indeks - 1]))[0] + '.fft')
            files.append(filepath)
            np.savetxt(filepath, tablica.view(float))
            indeks += 1

    def zwroc_tablice(self, file):
        """
        Metoda wypisująca do plików w postaci tablic, współczynniki Fouriera.
        Odpowiednio dla każdego pliku '*.png' w katalogu.
        :param lista_fft:
        :param path: Ścieżka
        :return: Pliki typu '*.txt' dla każdego obrazka '*.png'
        """
        fft = self.fft2d(ElementaryUnitCellImage().apply_for_desired_image(file))
        return fft.view(float)

    @staticmethod
    def fft1d(tablica):
        tab = np.fft.fftshift(np.fft.fft(tablica)) / len(tablica)
        return np.stack((tab.real, tab.imag), axis=-1)

    def wywolaj_fft1d(self, typ_struktury, repeat, len_num):
        """
        Metoda obliczająca wspoółczynniki Fouriera dla struktur jednowymiarowych.
        :param typ_struktury: Określa strukturę: periodyczną, Fibonacci, Thue-Morse'a
        :param repeat: Ilukrotnie powtórzyć kazdy element łańcucha. Potrzebne przy zwiększaniu ilości wekotorów sieci
        odwrotnej.
        :param len_num: Parametr, który decyduje z ilu segmentów składa się struktura. Pdaje się argument.
        :return: Współczynniki Fouriera w pliku tekstowym
        """
        # TODO: Sprawdzić normowanie
        if typ_struktury == 'TM':
            tab = ThueMorse(repeat, len_num).sequence()
            np.savetxt('tm_coef_' + str(repeat) + '*' + str(len_num) + '.fft', self.fft1d(tab))
        elif typ_struktury == 'F':
            tab = Fibonacci(repeat, len_num).sequence()
            np.savetxt('f_coef_' + str(repeat) + '*' + str(len_num) + '.fft', self.fft1d(tab))
        elif typ_struktury == 'P':
            tab = Periodic(repeat, len_num).sequence()
            np.savetxt('p_coef_' + str(repeat) + '*' + str(len_num) + '.fft', self.fft1d(tab))
        elif typ_struktury == 'C':
            tab = Heated(repeat).cos_sequence()
            np.savetxt('c_coef_' + str(repeat) + '.fft', self.fft1d(tab))
        elif typ_struktury == 'Custom':
            tab = Custom('real_Ni.fft').sequence()
            np.savetxt('heat_fft_Ni.fft', self.fft1d(tab))

        return np.fft.fftshift(np.fft.fft(tab))


if __name__ == "__main__":
    a = FFT().wywolaj_fft1d('P', 10, 377)

    #FFT().wypisz_do_pliku()
    #a = FFT().wypisz_do_pliku()
