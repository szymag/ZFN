# -*- coding: utf-8 -*-
from math import sqrt, exp, cosh
import numpy as np
from src.eig_problem.LoadFFT import LoadFFT2D
from src.eig_problem.ReciprocalVector import ReciprocalVector
from src.io.DataReader import ParsingData

from multiprocessing import Pool
import numexpr as ne


class EigenMatrix:
    class ReciprocalVectorGrid:
        def __init__(self, rec_vector_x, rec_vector_y):
            if rec_vector_y <= 0 or rec_vector_x <= 0:
                raise ValueError('Invalid arguments. Both should be positive.')
            self.rec_vector_x = rec_vector_x
            self.rec_vector_y = rec_vector_y

        def __str__(self):
            return str(self.rec_vector_x) + ' ' + str(self.rec_vector_x)

        def vectors_count(self):
            return self.rec_vector_y * self.rec_vector_x

        def coefficient_grid_size(self):
            return 2 * self.rec_vector_x - 1, 2 * self.rec_vector_y - 1

        def shift_to_middle_of_coeff_array(self):
            return np.array([self.rec_vector_x - 1, self.rec_vector_y - 1])

    def __init__(self, ReciprocalVectorGrid, vector_q,
                 input_parameters, material_A, material_B):
        # TODO: The constructor should be somehow modified to be more transparent
        # TODO: Info about materials shouldn't be here
        # TODO: **kwargs should allows to overwrite element in dictionary
        if isinstance(input_parameters, str):
            self.parameters = ParsingData(input_parameters)
        elif isinstance(input_parameters, dict):
            self.parameters = ParsingData(input_parameters)
        else:
            self.parameters = input_parameters

        self.material_A = self.parameters.material_constant(material_A)
        self.material_B = self.parameters.material_constant(material_B)
        self.ReciprocalVectorGrid = ReciprocalVectorGrid
        self.vectors_count = self.ReciprocalVectorGrid.vectors_count()
        self.gamma, self.mu0H0 = self.parameters.physical_constant()
        self.H0 = self.mu0H0 / self.parameters.mu0()
        self.tmp = LoadFFT2D(self.parameters.input_fft_file(),
                             self.ReciprocalVectorGrid.coefficient_grid_size())
        self.magnetization_sat = self.tmp.rescale_fourier_coefficient(self.material_A['Mo'], self.material_B['Mo'])
        self.exchange_len = self.tmp.rescale_fourier_coefficient(self.material_A['l'], self.material_B['l'])
        self.rec_vector_indexes = ReciprocalVector(self.vectors_count).lista_wektorow2d('min')
        self.shift_to_middle_of_coeff_array = ReciprocalVectorGrid.shift_to_middle_of_coeff_array()
        self.vector_q = vector_q

    def save_matrix_to_file(self):
        tmp = self.generate_and_fill_matrix()
        np.savetxt('matrix_to_eig_TheImpact_' + str(self.ReciprocalVectorGrid) +
                   '_q=[1e-9,0].txt', tmp.view(float))

    def generate_and_fill_matrix(self):
        matrix_to_eigen_problem = np.zeros((2 * self.vectors_count, 2 * self.vectors_count),
                                           dtype=complex)
        self.kroneker_delta(matrix_to_eigen_problem)
        w2 = self.rec_vector_indexes
        pool = Pool()
        tmp1 = np.array(pool.map(self.exchange_field, w2))
        tmp4 = np.array(pool.map(self.static_demagnetizing_field, w2))
        tmp2 = np.array(pool.map(self.dynamic_demagnetizing_field_in_plane, w2))
        tmp3 = np.array(pool.map(self.dynamic_demagnetizing_field_out_of_plane, w2))
        matrix_to_eigen_problem[self.vectors_count:, 0:self.vectors_count] += -tmp1 - tmp3 + tmp4
        matrix_to_eigen_problem[0:self.vectors_count, self.vectors_count:] += tmp1 + tmp2 - tmp4
        pool.close()
        return matrix_to_eigen_problem

    def kroneker_delta(self, matrix):
        for i in range(self.vectors_count, 2 * self.vectors_count):
            matrix[i - self.vectors_count][i] += 1.
            matrix[i][i - self.vectors_count] -= 1.

    def exchange_field(self, vector_2):
        div = self.parameters.lattice_const()[::-1]
        tab_from_wektor_1 = np.asfortranarray(
            np.broadcast_to(self.rec_vector_indexes, (self.vectors_count, self.vectors_count, 2)),
            dtype=int)
        tab_from_vec_l = np.transpose(tab_from_wektor_1, (1, 0, 2))

        vectors_diff = lambda a, b, s: ne.evaluate('{v1} - {v2} + {s}'.format(v1='a', v2='b', s='s'))

        tmp = vectors_diff(tab_from_wektor_1, tab_from_vec_l, self.shift_to_middle_of_coeff_array)
        tmp2 = vectors_diff(tab_from_vec_l, vector_2, self.shift_to_middle_of_coeff_array)

        tmp1 = self.magnetization_sat[tmp[:, :, 0], tmp[:, :, 1]]
        tmp3 = self.exchange_len[tmp2[:, :, 0], tmp2[:, :, 1]]

        vectors_sum = lambda a, b: ne.evaluate('2 * 3.14159265 * {v} + {q}'.format(v='a', q='b'))

        tmp4 = np.dot(
            vectors_sum(tab_from_vec_l / div, self.vector_q),
            vectors_sum(vector_2 / div, self.vector_q))

        e3 = lambda a, b, c, d: ne.evaluate('sum({a}*{b}*{c}/{d}, 0)'.format(a='a', b='b', c='c', d='d'))
        return e3(tmp1, tmp3, tmp4, self.H0)

    def dynamic_demagnetizing_field_in_plane(self, wektor_2):
        vec_2 = np.array(2 * np.pi * wektor_2) / self.parameters.lattice_const()[::-1]
        norm = sqrt((self.vector_q[0] + vec_2[0]) ** 2 + (self.vector_q[1] + vec_2[1]) ** 2)
        tmp1 = (self.vector_q[1] + vec_2[1]) ** 2 / norm ** 2
        tmp2 = 1 - cosh(norm * self.parameters.x()) * exp(-norm * self.parameters.thickness() / 2.)
        tmp3 = self.rec_vector_indexes - wektor_2 + self.shift_to_middle_of_coeff_array
        tmp4 = self.magnetization_sat[tmp3[:, 0], tmp3[:, 1]]
        e3 = lambda a, b, c, d: ne.evaluate('{a}*{b}*{c}/{d}'.format(a='a', b='b', c='c', d='d'))
        return e3(tmp1, tmp2, tmp4, self.H0)

    def dynamic_demagnetizing_field_out_of_plane(self, wektor_2):
        vec_2 = np.array(2 * np.pi * wektor_2) / self.parameters.lattice_const()[::-1]
        norm = sqrt((self.vector_q[0] + vec_2[0]) ** 2 + (self.vector_q[1] + vec_2[1]) ** 2)
        tmp1 = cosh(norm * self.parameters.x()) * exp(-norm * self.parameters.thickness() / 2.)
        tmp3 = self.rec_vector_indexes - wektor_2 + self.shift_to_middle_of_coeff_array
        tmp4 = self.magnetization_sat[tmp3[:, 0], tmp3[:, 1]]
        e3 = lambda a, b, c: ne.evaluate('{a}*{b}/{c}'.format(a='a', b='b', c='c'))
        return e3(tmp1, tmp4, self.H0)

    def static_demagnetizing_field(self, wektor_2):
        # along external field
        vector_1 = self.rec_vector_indexes
        vec_2 = np.array(2 * np.pi * wektor_2) / self.parameters.lattice_const()[::-1]
        vec_1 = np.array(2 * np.pi * vector_1) / self.parameters.lattice_const()[::-1]
        tmp1 = (vec_1[:, 0] - vec_2[0]) ** 2 / (np.linalg.norm(vec_1 - vec_2, axis=1) ** 2 + 1e-36)
        tmp2 = 1 - np.cosh(np.linalg.norm(vec_1 - vec_2, axis=1) * self.parameters.x()) * \
               np.exp(-np.linalg.norm(vec_1 - vec_2, axis=1) * self.parameters.thickness() / 2.)
        tmp3 = vector_1 - wektor_2 + self.shift_to_middle_of_coeff_array
        tmp4 = self.magnetization_sat[tmp3[:, 0], tmp3[:, 1]]
        e3 = lambda a, b, c, d: ne.evaluate('{a}*{b}*{c}/{d}'.format(a='a', b='b', c='c', d='d'))
        return e3(tmp1, tmp2, tmp4, self.H0)


if __name__ == "__main__":
    q = EigenMatrix(EigenMatrix.ReciprocalVectorGrid(3, 3), np.array([1e-9, 0]),
                    'InputParameter.yaml', 'Fe', 'Ni')
    q.save_matrix_to_file()
