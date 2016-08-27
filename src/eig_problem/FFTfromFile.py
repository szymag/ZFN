from math import sqrt
import pandas as pd
import numpy as np

class FFTfromFile:
    """
    Klasa, która wczytuje zadany plik tekstowy (domyślnie jest to 'fft1.txt' i wyciąga z niego informację o wektorach
    sieci odrwotnej wraz z odpowiadającymi im współczynnikami Fouriera.
    """

    def __init__(self, input_fft, tab_size):
        self.tmp_table = pd.read_csv(input_fft, delimiter=' ', dtype=float, header=None).values
        re = self.tmp_table[:, 0::2]
        im = self.tmp_table[:, 1::2] * 1j
        self.table = re + im
        self.tab_size = 2 * int(sqrt(tab_size))
        self.size = len(self.table[0])

    def coefficient(self):
        """
        Metoda tworząca listę współczynników. Na podstawie położenia w tablicy, określane jest położenie w liście
        :return: Lista współczynników.
        """
        return self.table[self.size // 2 - self.tab_size // 2:self.size // 2 + self.tab_size // 2 + 1,
               self.size // 2 - self.tab_size // 2:self.size // 2 + self.tab_size // 2 + 1]

    def fourier_coefficient(self, paramA, paramB):
        """
        Metoda ta tworzy słownik, gdzi kluczem jest wektor sieci odwrotnej, a wartością współczynnik Fouriera.
        :return: Słownik
        """
        k = self.tab_size // 2
        v = [(paramA - paramB) * i for i in self.coefficient()]
        v[k][k] += paramB
        assert v[k][k].imag == 0.
        return np.array(v)


if __name__ == "__main__":
    q = FFTfromFile('radius100.txt', 9)
    print(q.fourier_coefficient(1,2))
