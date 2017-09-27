from math import sqrt

import numpy as np
import pandas as pd

from src.eig_problem.WektorySieciOdwrotnej import WektorySieciOdwrotnej


class FFTfromFile2D:
    """
    Klasa, która wczytuje zadany plik tekstowy (domyślnie jest to 'fft1.txt' i wyciąga z niego informację o wektorach
    sieci odrwotnej wraz z odpowiadającymi im współczynnikami Fouriera.
    """

    def __init__(self):
        self.tmp_table = pd.read_csv(self.input_fft, delimiter=' ', dtype=float, header=None).values
        re = self.tmp_table[:, 0::2]
        im = self.tmp_table[:, 1::2] * 1j
        self.table = re + im
        self.coeff_number = self.ilosc_wektorow
        self.vector_max = WektorySieciOdwrotnej(self.a, self.b, self.coeff_number).lista_wektorow('max')

    def coefficient2d(self):
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
        v = [(self.MoA - self.MoB) * i for i in self.coefficient2d()]
        d = dict(zip(k, v))
        d[(0, 0)] = d[(0, 0)] + self.MoB
        assert d[(0, 0)].imag == 0.
        return d

    def exchange_length(self):
        """
        Metoda ta tworzy słownik, gdzi kluczem jest wektor sieci odwrotnej, a wartością współczynnik Fouriera.
        :return: Słownik
        """
        k = self.vector_max
        v = [(self.lA - self.lB) * i for i in self.coefficient2d()]
        d = dict(zip(k, v))
        d[(0, 0)] = d[(0, 0)] + self.lB
        return d


if __name__ == "__main__":
    q = FFTfromFile2D()
    print(q.fourier_coefficient()[(0, 0)])
