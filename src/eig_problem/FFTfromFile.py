from math import sqrt
import pandas as pd
import numpy as np

class FFTfromFile:
    """
    Klasa, która wczytuje zadany plik tekstowy (domyślnie jest to 'fft1.txt' i wyciąga z niego informację o wektorach
    sieci odrwotnej wraz z odpowiadającymi im współczynnikami Fouriera.
    """

    def __init__(self, input_fft, tab_size):

        """
        :type tab_size: tuple
        """
        self.tmp_table = pd.read_csv(input_fft, delimiter=' ', dtype=float, header=None).values
        re = self.tmp_table[:, 0::2]
        im = self.tmp_table[:, 1::2] * 1j
        self.table = re + im
        self.numbers_of_rec_vector = tab_size

    def coefficient(self):
        """
        Metoda tworząca listę współczynników. Na podstawie położenia w tablicy, określane jest położenie w liście
        :return: Lista współczynników.
        """
        middle_of_array = np.unravel_index(self.table.argmax(), self.table.shape)
        return self.table[ middle_of_array[0] - self.numbers_of_rec_vector[0] // 2:
               middle_of_array[0] + self.numbers_of_rec_vector[0] // 2 + self.numbers_of_rec_vector[0] % 2,
               middle_of_array[1] - self.numbers_of_rec_vector[1] // 2:
               middle_of_array[1] + self.numbers_of_rec_vector[1] // 2 + self.numbers_of_rec_vector[1] % 2]

    def fourier_coefficient(self, paramA, paramB):
        """
        Metoda ta tworzy słownik, gdzi kluczem jest wektor sieci odwrotnej, a wartością współczynnik Fouriera.
        :return: Słownik
        """
        k = self.numbers_of_rec_vector // 2
        v = [(paramA - paramB) * i for i in self.coefficient()]
        v[k][k] += paramB
        assert v[k][k].imag == 0.
        return np.array(v)


if __name__ == "__main__":
    q = FFTfromFile('penrose.txt', (4,4))
    print(q.coefficient())
