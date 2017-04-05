import numpy as np

from src.eig_problem.ParametryMaterialowe import ParametryMaterialowe


class FFTfromFile1D:
    def __init__(self, input_fft):
        self.file = np.transpose(np.loadtxt(input_fft))
        self.coeff = self.file[0] + self.file[1] * 1j
        self.ilosc_wektorow = len(self.file[0]) // 2

    def coefficient1d(self):
        index = len(self.file[0]) // 2
        return np.array(self.coeff[index - int(self.ilosc_wektorow) + 1:index + int(self.ilosc_wektorow)])

    def fourier_coefficient(self, paramA, paramB):
        """
        Metoda ta tworzy słownik, gdzi kluczem jest wektor sieci odwrotnej, a wartością współczynnik Fouriera.
        :return: Słownik
        """
        index0 = self.ilosc_wektorow - 1
        tab = (paramA - paramB) * self.coefficient1d()
        tab[index0] += paramB
        assert tab[index0].imag == 0.
        return tab

