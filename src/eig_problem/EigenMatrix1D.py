# -*- coding: utf-8 -*-
import numpy as np
from math import radians

from src.eig_problem.LoadFFT import LoadFFT1D
from src.eig_problem.ReciprocalVector import ReciprocalVector
from src.io.DataReader import ParsingData


class EigenMatrix1D:
    def __init__(self, bloch_vec, input_parameters, material_A, material_B):
        if isinstance(input_parameters, str):
            self.parameters = ParsingData(input_parameters)
        elif isinstance(input_parameters, dict):
            self.parameters = ParsingData(input_parameters)
        else:
            self.parameters = input_parameters

        self.material_A = self.parameters.material_constant(material_A)
        self.material_B = self.parameters.material_constant(material_B)

        self.gamma, self.mu0H0 = self.parameters.physical_constant()
        self.H0 = self.mu0H0 / self.parameters.mu0()
        self.tmp = LoadFFT1D(self.parameters.input_fft_file())
        self.vectors_count = self.tmp.vectors_count
        self.magnetization_sat = self.tmp.fourier_coefficient(self.material_A['Mo'], self.material_B['Mo'])
        self.exchange_len = self.tmp.fourier_coefficient(self.material_A['l'], self.material_B['l'])
        self.reciprocal_vec = ReciprocalVector(self.vectors_count).lista_wektorow1d('min')
        self.shift_to_middle_of_coeff_array = len(self.reciprocal_vec) - 1
        self.bloch_vec = bloch_vec

    def generate_and_fill_matrix(self):
        matrix = np.zeros((2 * self.vectors_count, 2 * self.vectors_count), dtype=complex)
        indeks = self.vectors_count
        self.kroneker_delta(matrix)
        for i in range(indeks, 2 * indeks):
            w1 = self.reciprocal_vec[i - indeks]
            w2 = self.reciprocal_vec
            ex = self.exchange_field(w1, w2)
            dyn_in_plane = self.dynamic_demagnetizing_field_in_plane(w1, w2)
            dyn_out_plane = self.dynamic_demagnetizing_field_out_of_plane(w1, w2)
            static = self.static_demagnetizing_field(w1, w2)
            matrix[i][np.arange(indeks)] += -ex - dyn_in_plane * np.sin\
                (radians(self.parameters.angle())) ** 2 + static * np.cos(radians(self.parameters.angle())) ** 2  # yx
            matrix[i - indeks][np.arange(indeks, 2 * indeks)] += \
                ex + dyn_out_plane - static * np.cos(radians(self.parameters.angle()))**2  # xy
        return matrix

    def exp_function(self, vec_1, vec_2):
        return np.exp(-abs((vec_1 + vec_2)) * self.parameters.thickness() / 2)

    def kroneker_delta(self, matrix):
        for i in range(self.vectors_count, 2 * self.vectors_count):
            matrix[i - self.vectors_count][i] += 1.
            matrix[i][i - self.vectors_count] -= 1.

    def exchange_field(self, vec_1, vec_2):
        vec_l = np.transpose(np.broadcast_to(self.reciprocal_vec, (self.vectors_count, self.vectors_count)))
        tmp1 = self.exchange_len[vec_l - vec_2 + self.shift_to_middle_of_coeff_array]
        tmp2 = self.magnetization_sat[vec_1 - vec_l + self.shift_to_middle_of_coeff_array]
        tmp3 = (self.bloch_vec + 2 * np.pi * vec_2 /
                self.parameters.lattice_const()[0]) * (2 * np.pi * vec_l /
                                                       self.parameters.lattice_const()[0] + self.bloch_vec)
        return np.sum(tmp1 * tmp2 * tmp3, axis=0) / self.H0

    def dynamic_demagnetizing_field_in_plane(self, vec_1, vec_2):
        tmp3 = self.magnetization_sat[vec_1 - vec_2 + self.shift_to_middle_of_coeff_array]
        tmp2 = 1 - self.exp_function(self.bloch_vec, (2 * np.pi * vec_2 / self.parameters.lattice_const()[0]))
        return tmp3 * tmp2 / self.H0

    def dynamic_demagnetizing_field_out_of_plane(self, vec_1, vec_2):
        tmp3 = self.magnetization_sat[vec_1 - vec_2 + self.shift_to_middle_of_coeff_array]
        tmp2 = self.exp_function(self.bloch_vec, (2 * np.pi * vec_2 / self.parameters.lattice_const()[0]))
        return tmp2 * tmp3 / self.H0

    def static_demagnetizing_field(self, vec_1, vec_2):
        tmp1 = self.magnetization_sat[vec_1 - vec_2 + self.shift_to_middle_of_coeff_array]
        co = np.cosh(abs((2 * np.pi * vec_1 / self.parameters.lattice_const()[0] - 2 * np.pi * vec_2 /
                          self.parameters.lattice_const()[0])) * self.parameters.x())
        tmp2 = 1 - co*np.exp(-abs((2 * np.pi * vec_1 / self.parameters.lattice_const()[0] - 2 * np.pi * vec_2 /
                                   self.parameters.lattice_const()[0])) * self.parameters.thickness() / 2)
        tmp1[self.shift_to_middle_of_coeff_array] = 0
        return tmp2 * tmp1 / self.H0

    def save_matrix_to_file(self):
        self.generate_and_fill_matrix(1e-9)
        np.savetxt('macierz.txt', np.array(self.matrix))


if __name__ == "__main__":
    q = EigenMatrix1D('soko.fft', 1e-9)
    q.save_matrix_to_file()
