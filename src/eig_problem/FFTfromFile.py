from math import sqrt

import numpy as np

from src.eig_problem.ParametryMaterialowe import ParametryMaterialowe
from src.eig_problem.WektorySieciOdwrotnej import WektorySieciOdwrotnej


class FFTfromFile(ParametryMaterialowe):
    """
    Klasa, która wczytuje zadany plik tekstowy (domyślnie jest to 'fft1.txt' i wyciąga z niego informację o wektorach
    sieci odrwotnej wraz z odpowiadającymi im współczynnikami Fouriera.
    """

    def __init__(self, ilosc_wektorow, typ_pole_wymiany, filepath='fft1.txt'):
        ParametryMaterialowe.__init__(self, ilosc_wektorow, typ_pole_wymiany)
        self.table = np.loadtxt(filepath).view(complex)
        self.coeff_number = ilosc_wektorow
        self.vector_max = WektorySieciOdwrotnej(self.a, self.b, self.coeff_number).lista_wektorow('max')

    def coefficient(self):
        # TODO: Poprawić wybieranie współczynników
        """
        Metoda tworząca listę współczynników. Na podstawie położenia w tablicy, określane jest położenie w liście
        :return: Lista współczynników.
        """
        index = (sqrt(len(self.vector_max)) - 1.) / 2.
        index1 = int(len(self.table) / 2.) - index
        index2 = int(sqrt(len(self.vector_max)))
        coefficient = list(np.zeros(index2 * index2))
        for i in range(index2):
            for j in range(index2):
                coefficient[j + i * index2] = self.table[i + index1][j + index1]
        return coefficient

    def fourier_coefficient(self):
        """
        Metoda ta tworzy słownik, gdzi kluczem jest wektor sieci odwrotnej, a wartością współczynnik Fouriera.
        :return: Słownik
        """
        k = self.vector_max
        v = [(self.MoCo - self.MoPy) * i for i in self.coefficient()]
        d = dict(zip(k, v))
        d[(0, 0)] = d[(0, 0)] + self.MoPy
        assert d[(0, 0)].imag == 0.
        return d

    def exchange_length(self):
        """
        Metoda ta tworzy słownik, gdzi kluczem jest wektor sieci odwrotnej, a wartością współczynnik Fouriera.
        :return: Słownik
        """
        k = self.vector_max
        v = [(self.lCo - self.lPy) * i for i in self.coefficient()]
        d = dict(zip(k, v))
        d[(0, 0)] = d[(0, 0)] + self.lPy
        return d


if __name__ == "__main__":
    q = FFTfromFile(121, 'I')
    print(q.fourier_coefficient()[(0, 0)])
