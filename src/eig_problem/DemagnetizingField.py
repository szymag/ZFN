import numpy as np
from math import cos
from src.eig_problem.ParametryMaterialowe import ParametryMaterialowe
from src.eig_problem.WektorySieciOdwrotnej import WektorySieciOdwrotnej
from src.eig_problem.FFTfromFile1D import FFTfromFile1D
import matplotlib.pyplot as plt


class StaticDemagnetizingField_1D:
    def __init__(self, input_fft, field_angle):
        self.tmp = FFTfromFile1D(input_fft)
        self.field_angle = field_angle
        self.ilosc_wektorow = self.tmp.ilosc_wektorow
        self.input_fft = input_fft
        self.coefficient = self.tmp.fourier_coefficient(ParametryMaterialowe.MoA, ParametryMaterialowe.MoB)
        self.reciprocal_vector = 2 * np.pi * WektorySieciOdwrotnej(self.ilosc_wektorow).lista_wektorow1d('max') \
                                 / ParametryMaterialowe.a

    def demagnetizing_field_at_point(self, location):
        tmp1 = np.ones(len(self.coefficient)) - np.exp(-np.absolute(self.reciprocal_vector) * ParametryMaterialowe.d / 2.)
        return np.sum(tmp1 * self.coefficient * np.exp(1j * location * self.reciprocal_vector)).real \
               * cos(self.field_angle / 360 * 2 * np.pi)

    def demagnetizing_field(self):
        elementary_cell = np.linspace(-ParametryMaterialowe.a, ParametryMaterialowe.a, 1000)
        demag = np.zeros(len(elementary_cell))
        for i in enumerate(elementary_cell):
            demag[i[0]] = self.demagnetizing_field_at_point(i[1])
        return elementary_cell, demag

    def demagnetizing_field_plot(self):
        tmp = self.demagnetizing_field()
        plt.plot(tmp[0], tmp[1])
        plt.ylim([-120000, 120000])
        plt.show()

if __name__ == "__main__":
    q = StaticDemagnetizingField_1D('c_coef_1000.txt', 0).demagnetizing_field_plot()
