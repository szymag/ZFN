import pandas as pd
import numpy as np
import sys

class LoadFFT:
    def __init__(self, input_fft, tab_size):
        try:
            self.loaded_fourier_coefficient_from_file = \
                pd.read_csv(input_fft, delimiter=' ', dtype=float, header=None).values
        except IOError:
            print('Can not find such file. Generate fourier coefficient in FFT class in fft_from_file directory')
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


if __name__ == "__main__":
    pass
