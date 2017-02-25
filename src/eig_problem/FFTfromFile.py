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


    def fourier_coefficient(self, paramA, paramB):
        v = (paramA - paramB) * self.coefficient()
        middle_of_array = self.middle_of_array(v)
        v[middle_of_array[0]][middle_of_array[1]] += paramB
        assert v[middle_of_array[0]][middle_of_array[1]].imag == 0.
        return np.array(v)

    def coefficient(self):
        middle_of_array = self.middle_of_array(np.array(self.table))
        return self.table[middle_of_array[0] - self.numbers_of_rec_vector[0] // 2:
               middle_of_array[0] + self.numbers_of_rec_vector[0] // 2 + self.numbers_of_rec_vector[0] % 2,
               middle_of_array[1] - self.numbers_of_rec_vector[1] // 2:
               middle_of_array[1] + self.numbers_of_rec_vector[1] // 2 + self.numbers_of_rec_vector[1] % 2]

    def middle_of_array(self, array):
        return np.unravel_index(abs(array).argmax(), array.shape)


if __name__ == "__main__":
    pass
