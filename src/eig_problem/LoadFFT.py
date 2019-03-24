import pandas as pd
import numpy as np
import sys


class LoadFFT2D:
    def __init__(self, input_fft, tab_size):
        try:
            self.loaded_fourier_coefficient_from_file = \
                pd.read_csv(input_fft, delimiter=' ', dtype=float, header=None).values
        except IOError:
            print('Can not find such file: ' + input_fft + ' \n Generate fourier coefficient in FFT class in fft_from_file directory')
            sys.exit()
        self.tab_size = tab_size

    def rescale_fourier_coefficient(self, paramA, paramB):
        v = (paramA - paramB) * self.crop_table()
        middle_of_array = self.return_middle_of_array(v)
        v[middle_of_array[0]][middle_of_array[1]] += paramB
        assert v[middle_of_array[0]][middle_of_array[1]].imag == 0.
        return np.array(v)

    def crop_table(self):
        table_to_crop = self.fourier_coefficient()
        middle_of_array = self.return_middle_of_array(table_to_crop)
        return table_to_crop[middle_of_array[0] - self.tab_size[0] // 2:
               middle_of_array[0] + self.tab_size[0] // 2 + self.tab_size[0] % 2,
               middle_of_array[1] - self.tab_size[1] // 2:
               middle_of_array[1] + self.tab_size[1] // 2 + self.tab_size[1] % 2]

    def return_middle_of_array(self, array):
        return np.array(np.unravel_index(abs(array).argmax(), array.shape))

    def fourier_coefficient(self):
        re = self.loaded_fourier_coefficient_from_file[:, 0::2]
        im = self.loaded_fourier_coefficient_from_file[:, 1::2] * 1j
        table = np.array(re + im, dtype=complex)
        return table


class LoadFFT1D:
    def __init__(self, input_fft):
        try:
            self.file = np.transpose(np.loadtxt(input_fft))
        except IOError:
            print('Can not find such file. Generate fourier coefficient in FFT class in fft_from_file directory')
            sys.exit()
        self.coeff = self.file[0] + self.file[1] * 1j
        self.vectors_count = len(self.file[0]) // 2

    def coefficient1d(self):
        index = len(self.file[0]) // 2
        return np.array(self.coeff[index - int(self.vectors_count) + 1:index + int(self.vectors_count)])

    def fourier_coefficient(self, paramA, paramB):
        """
        Metoda ta tworzy słownik, gdzi kluczem jest wektor sieci odwrotnej, a wartością współczynnik Fouriera.
        :return: Słownik
        """
        index0 = self.vectors_count - 1
        tab = (paramA - paramB) * self.coefficient1d()
        tab[index0] += paramB
        assert tab[index0].imag == 0.
        return tab


if __name__ == "__main__":
    pass
