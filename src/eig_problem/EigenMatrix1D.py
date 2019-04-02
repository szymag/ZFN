# -*- coding: utf-8 -*-
import numpy as np
from math import radians, sin
from collections.abc import Iterable
from multiprocessing import Pool
import numexpr as ne

from src.eig_problem.LoadFFT import LoadFFT1D
from src.eig_problem.ReciprocalVector import ReciprocalVector
from src.io.DataReader import ParsingData
from src.utils.cProfiler import do_cprofile


class EigenMatrix1D:
    def __init__(self, bloch_vec, input_parameters, material_A, material_B, bloch_vec_perp=0):
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
        self.tmp = LoadFFT1D(self.parameters.fft_data())
        self.vectors_count = self.tmp.vectors_count
        self.magnetization_sat = self.tmp.fourier_coefficient(self.material_A['Mo'], self.material_B['Mo'])
        self.exchange_len = self.tmp.fourier_coefficient(self.material_A['l'], self.material_B['l'])
        self.reciprocal_vec = ReciprocalVector(self.vectors_count).lista_wektorow1d('min')
        self.shift_to_middle_of_coeff_array = len(self.reciprocal_vec) - 1
        if isinstance(bloch_vec, Iterable):
            self.bloch_vec = bloch_vec[0]
            self.bloch_vec_perp = bloch_vec[1]
        else:
            self.bloch_vec = bloch_vec
            self.bloch_vec_perp = bloch_vec_perp

    @do_cprofile
    def generate_and_fill_matrix(self):
        matrix = np.zeros((2 * self.vectors_count, 2 * self.vectors_count), dtype=np.complex128)
        indeks = self.vectors_count
        self.kroneker_delta(matrix)

        w1 = self.reciprocal_vec
        pool = Pool()
        ex = np.array(pool.map(self.exchange_field, w1))
        dyn_in_plane = np.array(pool.map(self.dynamic_demagnetizing_field_in_plane, w1))
        dyn_out_plane = np.array(pool.map(self.dynamic_demagnetizing_field_out_of_plane, w1))
        static =  np.array(pool.map(self.static_demagnetizing_field, w1))
        matrix[indeks:, 0:indeks] += -ex - dyn_in_plane * np.sin(radians(self.parameters.angle())) ** 2 + static * np.cos(radians(self.parameters.angle())) ** 2  # yx
        matrix[0:indeks, indeks:] += ex + dyn_out_plane - static * np.cos(radians(self.parameters.angle()))**2  # xy
        pool.close()

        return matrix

    def exp_function(self, vec_1, vec_2):
        tmp = self.demag_field_factor(vec_1, vec_2) # if disabled, spin wave spectra are calculated on surface
        return tmp * np.exp(-np.linalg.norm([vec_1 + vec_2, self.bloch_vec_perp]) * self.parameters.thickness() / 2)

    def demag_field_factor(self, vec_1, vec_2):
        """
        This factor is says where spin wave spectra are calculated. By default it is done on
        surface so we can neglect this term : x=0
        """
        x = np.linspace(-self.parameters.thickness()/2, self.parameters.thickness()/2, 50)
        norm_vec = np.cosh(np.linalg.norm([vec_1 + vec_2, self.bloch_vec_perp]) * x[:, np.newaxis])
        tmp = np.trapz(norm_vec, x, axis=0) / self.parameters.thickness()
        #print(tmp*sin(radians(self.parameters.perpendicular_bloch_vector()[0]))**2 / 2)
        return tmp

    def kroneker_delta(self, matrix):
        for i in range(self.vectors_count, 2 * self.vectors_count):
            matrix[i - self.vectors_count][i] += 1.
            matrix[i][i - self.vectors_count] -= 1.

    def exchange_field(self, vec_1):
        mid_point = self.shift_to_middle_of_coeff_array
        rec_vec = self.reciprocal_vec
        lat_const = self.parameters.lattice_const()[0]
        bloch_vec = self.bloch_vec
        bloch_vec_perp = self.bloch_vec_perp
        pi = np.pi
        H0 = self.H0
        vec_l = np.broadcast_to(rec_vec, (self.vectors_count, self.vectors_count)).T
        tmp1 = self.exchange_len[ne.evaluate('vec_l - rec_vec + mid_point')]
        tmp2 = self.magnetization_sat[ne.evaluate('vec_1 - vec_l + mid_point')]
        # tmp3 = ne.evaluate('(bloch_vec + 2 * pi * rec_vec / lat_const) * (2 * pi * vec_l / lat_const + bloch_vec) + bloch_vec_perp**2')

        # e3 = lambda a, b, c, d: ne.evaluate('sum({a}*{b}*{c}/{d}, 0)'.format(a='a', b='b', c='c', d='d'))
        # return e3(tmp1, tmp2, tmp3, self.H0)

        return ne.evaluate('sum(tmp1*tmp2*((bloch_vec + 2 * pi * rec_vec / lat_const) * (2 * pi * vec_l / lat_const + bloch_vec) + bloch_vec_perp**2)/H0, 0)')

    def exchange_field_simplified(self, vec_1):
        tmp1 = self.exchange_len[vec_1 - self.reciprocal_vec + self.shift_to_middle_of_coeff_array]
        tmp3 = (self.bloch_vec + 2 * np.pi * self.reciprocal_vec /
                self.parameters.lattice_const()[0]) * (2 * np.pi * vec_1 /
                                                       self.parameters.lattice_const()[0] + self.bloch_vec)
        return tmp1 * tmp3 / self.H0

    def dynamic_demagnetizing_field_in_plane(self, vec_1):
        tmp3 = self.magnetization_sat[vec_1 - self.reciprocal_vec + self.shift_to_middle_of_coeff_array]
        tmp2 = 1 - self.exp_function(self.bloch_vec,
                                     (2 * np.pi * self.reciprocal_vec / self.parameters.lattice_const()[0]))
        tmp1 = (self.reciprocal_vec + self.bloch_vec)**2 / np.linalg.norm([self.reciprocal_vec + self.bloch_vec, self.bloch_vec_perp])**2
        e3 = lambda a, b, c, d: ne.evaluate('{a}*{b}*{c}/{d}'.format(a='a', b='b', c='c', d='d'))
        return e3(tmp3, tmp2, tmp1, self.H0)
        #return tmp3 * tmp2 * tmp1 / self.H0

    def dynamic_demagnetizing_field_out_of_plane(self, vec_1):
        tmp3 = self.magnetization_sat[vec_1 - self.reciprocal_vec + self.shift_to_middle_of_coeff_array]
        tmp2 = self.exp_function(self.bloch_vec, (2 * np.pi * self.reciprocal_vec / self.parameters.lattice_const()[0]))
        e3 = lambda a, b, c: ne.evaluate('{a}*{b}/{c}'.format(a='a', b='b', c='c'))
        return e3(tmp2, tmp3, self.H0)
        #return tmp2 * tmp3 / self.H0

    def static_demagnetizing_field(self, vec_1):
        tmp1 = self.magnetization_sat[vec_1 - self.reciprocal_vec + self.shift_to_middle_of_coeff_array]
        co = np.cosh(abs((2 * np.pi * vec_1 / self.parameters.lattice_const()[0] - 2 * np.pi * self.reciprocal_vec /
                          self.parameters.lattice_const()[0])) * self.parameters.x())
        tmp2 = 1 - co*np.exp(-abs((2 * np.pi * vec_1 / self.parameters.lattice_const()[0] - 2 * np.pi * self.reciprocal_vec /
                                   self.parameters.lattice_const()[0])) * self.parameters.thickness() / 2)
        e3 = lambda a, b, c: ne.evaluate('{a}*{b}/{c}'.format(a='a', b='b', c='c'))
        return e3(tmp2, tmp1, self.H0)
        #return tmp2 * tmp1 / self.H0

    def save_matrix_to_file(self):
        np.savetxt('macierz.txt', np.array(self.generate_and_fill_matrix()).view(float))


if __name__ == "__main__":
    a = EigenMatrix1D([0, 1e-9], './test.yaml', 'Co', 'Py')
    a.generate_and_fill_matrix()
    #a.save_matrix_to_file()
