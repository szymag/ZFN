import numpy as np
from math import cos
from src.eig_problem.ParametryMaterialowe import ParametryMaterialowe
from src.eig_problem.WektorySieciOdwrotnej import WektorySieciOdwrotnej
from src.eig_problem.FFTfromFile1D import FFTfromFile1D
import matplotlib.pyplot as plt


class StaticDemagnetizingField1D:
    def __init__(self, input_fft, field_angle):
        self.tmp = FFTfromFile1D(input_fft)
        self.field_angle = field_angle
        self.ilosc_wektorow = self.tmp.ilosc_wektorow
        self.input_fft = input_fft
        self.coefficient = self.tmp.fourier_coefficient(ParametryMaterialowe.MoA, ParametryMaterialowe.MoB)
        self.reciprocal_vector = 2 * np.pi * WektorySieciOdwrotnej(self.ilosc_wektorow).lista_wektorow1d('max') \
                                 / ParametryMaterialowe.a

    def elementary_cell_reconstruction(self, grid):
        coefficient = np.transpose(np.loadtxt('c_coef_100.txt').view(complex))
        x = np.linspace(-ParametryMaterialowe.a, 0, grid)
        tmp = np.zeros(grid)
        for ind in enumerate(x):
            tmp[ind[0]] = abs(self.inverse_discrete_fourier_transform(coefficient, ind[1]))
        return x , tmp / (10 / 7) + 0.3

    def inverse_discrete_fourier_transform(self, data, vector_position):
        reciprocal_vectors = np.array(2 * np.pi * WektorySieciOdwrotnej(max(data.shape)).lista_wektorow1d('min')
                                      / ParametryMaterialowe.a)
        return np.sum(data * np.exp(1j * reciprocal_vectors * vector_position))

    def demagnetizing_field_at_point(self, location):
        tmp1 = np.ones(len(self.coefficient)) - np.exp(-np.absolute(self.reciprocal_vector) * ParametryMaterialowe.d / 2.)
        return np.sum(tmp1 * self.coefficient * np.exp(1j * location * self.reciprocal_vector)).real \
               * cos(self.field_angle / 360 * 2 * np.pi)

    def demagnetizing_field(self):
        elementary_cell = np.linspace(-ParametryMaterialowe.a, 0, 1000)
        demag = np.zeros(len(elementary_cell))
        for i in enumerate(elementary_cell):
            demag[i[0]] = self.demagnetizing_field_at_point(i[1])
        return elementary_cell, demag / np.max(demag)

    def demagnetizing_field_plot(self):
        tmp = self.demagnetizing_field()
        tmp1 = self.elementary_cell_reconstruction(500)
        plt.plot(tmp1[0], tmp1[1])
        plt.plot(tmp[0], tmp[1])
        #plt.ylim([-120000, 120000])
        plt.savefig('normalized_demagnetized_field_sin_like.png')
        plt.show()

if __name__ == "__main__":
    q = StaticDemagnetizingField1D('c_coef_100.txt', 0).demagnetizing_field_plot()
