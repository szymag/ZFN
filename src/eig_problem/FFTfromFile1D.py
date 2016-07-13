import numpy as np
from src.eig_problem.ParametryMaterialowe import ParametryMaterialowe
from src.eig_problem.WektorySieciOdwrotnej import WektorySieciOdwrotnej
import matplotlib.pyplot as plt

class FFTfromFile1D(ParametryMaterialowe):
    def __init__(self):
        ParametryMaterialowe.__init__(self)
        self.file = np.transpose(np.loadtxt(self.input_fft))
        self.coeff = self.file[0] + self.file[1] * 1j
        self.vector_max = WektorySieciOdwrotnej(self.a, self.b, self.ilosc_wektorow).lista_wektorow1d('max')

    def coefficient1d(self):
        index = len(self.coeff) // 2
        return self.coeff[index  - int(self.ilosc_wektorow) + 1:index + int(self.ilosc_wektorow)]

    def fourier_coefficient(self):
        """
        Metoda ta tworzy słownik, gdzi kluczem jest wektor sieci odwrotnej, a wartością współczynnik Fouriera.
        :return: Słownik
        """
        k = self.vector_max
        v = [(self.MoA - self.MoB) * i for i in self.coefficient1d()]
        d = dict(zip(k, v))
        d[0] = d[0] + self.MoB
        assert d[0].imag == 0.
        return d

    def exchange_length(self):
        """
        Metoda ta tworzy słownik, gdzi kluczem jest wektor sieci odwrotnej, a wartością współczynnik Fouriera.
        :return: Słownik
        """
        k = self.vector_max
        v = [(self.lA - self.lB) * i for i in self.coefficient1d()]
        d = dict(zip(k, v))
        d[0] = d[0] + self.lB
        return d

    def reconstruction(self, r):
        coeff = self.coefficient1d()
        ex = np.exp(1j * self.vector_max * r / self.a)
        return abs(np.sum(coeff * ex))

    def plot(self):
        tab = np.zeros(len(self.vector_max))
        print(len(self.vector_max))
        plt.plot(np.arange(0, len(self.vector_max)), abs(self.coefficient1d()))
        plt.show()


