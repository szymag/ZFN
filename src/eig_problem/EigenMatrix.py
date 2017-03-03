# -*- coding: utf-8 -*-
from math import sqrt, exp, cosh
import numpy as np
from src.eig_problem.FFTfromFile import FFTfromFile
from src.eig_problem.InputParameter import InputParameter
from src.eig_problem.ReciprocalVector import ReciprocalVector
from multiprocessing import Pool
import numexpr as ne


class EigenMatrix:
    def __init__(self, input_fft, rec_vector_x, rec_vector_y, wektor_q, a=InputParameter.a, b=InputParameter.b,
                 MoA=InputParameter.MoA, MoB=InputParameter.MoB, lA=InputParameter.lA,
                 lB=InputParameter.lB, d=InputParameter.d,
                 x=InputParameter.x, H0=InputParameter.H0):
        self.lattice_const_x = a
        self.lattice_const_y = b
        self.thickness = d
        self.rec_vector_x = rec_vector_x
        self.rec_vector_y = rec_vector_y
        self.amount_rec_vector = rec_vector_x * rec_vector_y
        self.x = x
        self.H0 = H0
        self.tmp = FFTfromFile(input_fft, (2*self.rec_vector_x-1, 2*self.rec_vector_y-1))
        self.magnetization_sat = self.tmp.fourier_coefficient(MoA, MoB)
        self.exchange_len = self.tmp.fourier_coefficient(lA, lB)
        self.rec_vector_indexes = ReciprocalVector(self.amount_rec_vector).lista_wektorow2d('min')
        #
        self.shift_to_middle_of_coeff_array = np.array([self.rec_vector_x - 1, self.rec_vector_y - 1])
        self.wektor_q = wektor_q

    def save_matrix_to_file(self):
        tmp = self.generate_and_fill_matrix()
        np.savetxt('matrix_to_eig_TheImpact_' + str(self.rec_vector_x) + ' ' + str(self.rec_vector_x) +
                   '_q=[1e-9,0].txt', tmp.view(float))

    def generate_and_fill_matrix(self):
        self.macierz_M = np.zeros((2 * self.amount_rec_vector, 2 * self.amount_rec_vector),
                                  dtype=complex)
        self.kroneker_delta(self.macierz_M)
        w2 = self.rec_vector_indexes
        pool = Pool()
        tmp1 = np.array(pool.map(self.exchange_field, w2))
        tmp4 = np.array(pool.map(self.static_demagnetizing_field, w2))
        tmp2 = np.array(pool.map(self.dynamic_demagnetizing_field_in_plane, w2))
        tmp3 = np.array(pool.map(self.dynamic_demagnetizing_field_out_of_plane, w2))
        self.macierz_M[self.amount_rec_vector:, 0:self.amount_rec_vector] += -tmp1 - tmp3 + tmp4
        self.macierz_M[0:self.amount_rec_vector, self.amount_rec_vector:] += tmp1 + tmp2 - tmp4
        pool.close()
        return self.macierz_M

    def kroneker_delta(self, matrix):
        for i in range(self.amount_rec_vector, 2 * self.amount_rec_vector):
            matrix[i - self.amount_rec_vector][i] += 1.
            matrix[i][i - self.amount_rec_vector] -= 1.

    # noinspection PyTypeChecker
    def exchange_field(self, wektor_2):
        div = [self.lattice_const_y, self.lattice_const_x]
        H0 = self.H0
        tab_from_wektor_1 = np.asfortranarray(
            np.broadcast_to(self.rec_vector_indexes, (self.amount_rec_vector , self.amount_rec_vector, 2)),
            dtype=int)
        tab_from_vec_l = np.transpose(tab_from_wektor_1, (1, 0, 2))

        e1 = lambda a, b, s: ne.evaluate(
                '{v1} - {v2} + {s}'.format(v1='a', v2='b', s='s'))

        tmp = e1(tab_from_wektor_1, tab_from_vec_l, self.shift_to_middle_of_coeff_array)
        tmp2 = e1(tab_from_vec_l, wektor_2, self.shift_to_middle_of_coeff_array)

        tmp1 = self.magnetization_sat[tmp[:, :, 0], tmp[:, :, 1]]
        tmp3 = self.exchange_len[tmp2[:, :, 0], tmp2[:, :, 1]]

        e2 = lambda a, b: ne.evaluate(
                '2 * 3.14159265 * {v} + {q}'.format(v='a', q='b'))

        tmp4 = np.dot(
                e2(tab_from_vec_l / div, self.wektor_q),
                e2(wektor_2 / div, self.wektor_q))

        return ne.evaluate('sum(tmp1*tmp3*tmp4/H0, 0)')

    # noinspection PyTypeChecker
    def dynamic_demagnetizing_field_in_plane(self, wektor_2):
        H0 = self.H0
        wekt_2 = np.array(2 * np.pi * wektor_2 / [self.lattice_const_y, self.lattice_const_x])
        norm = sqrt((self.wektor_q[0] + wekt_2[0]) ** 2 + (self.wektor_q[1] + wekt_2[1]) ** 2)
        tmp1 = (self.wektor_q[1] + wekt_2[1]) ** 2 / norm ** 2
        tmp2 = 1 - cosh(norm * self.x) * exp(-norm * self.thickness / 2.)
        tmp3 = self.rec_vector_indexes - wektor_2 + self.shift_to_middle_of_coeff_array
        tmp4 = self.magnetization_sat[tmp3[:, 0], tmp3[:, 1]]
        return ne.evaluate('tmp1 * tmp2 * tmp4 / H0')

    # noinspection PyTypeChecker
    def dynamic_demagnetizing_field_out_of_plane(self, wektor_2):
        H0 = self.H0
        wekt_2 = np.array(2 * np.pi * wektor_2 / [self.lattice_const_y, self.lattice_const_x])
        norm = sqrt((self.wektor_q[0] + wekt_2[0]) ** 2 + (self.wektor_q[1] + wekt_2[1]) ** 2)
        tmp1 = cosh(norm * self.x) * exp(-norm * self.thickness / 2.)
        tmp3 = self.rec_vector_indexes - wektor_2 + self.shift_to_middle_of_coeff_array
        tmp4 = self.magnetization_sat[tmp3[:, 0], tmp3[:, 1]]
        return ne.evaluate('tmp1 * tmp4 / H0')

    # noinspection PyTypeChecker
    def static_demagnetizing_field(self, wektor_2):
        # along external field
        H0 = self.H0
        wektor_1 = self.rec_vector_indexes
        wekt_2 = np.array(2 * np.pi * wektor_2 / [self.lattice_const_y, self.lattice_const_x])
        wekt_1 = 2 * np.pi * wektor_1 / [self.lattice_const_y, self.lattice_const_x]
        tmp1 = (wekt_1[:, 0] - wekt_2[0]) ** 2 / (np.linalg.norm(wekt_1 - wekt_2, axis=1) ** 2 + 1e-36)
        tmp2 = 1 - np.cosh(np.linalg.norm(wekt_1 - wekt_2, axis=1) * self.x) * \
                   np.exp(-np.linalg.norm(wekt_1 - wekt_2, axis=1) * self.thickness / 2.)
        tmp3 = wektor_1 - wektor_2 + self.shift_to_middle_of_coeff_array
        tmp4 = self.magnetization_sat[tmp3[:, 0], tmp3[:, 1]]
        return ne.evaluate('tmp1 * tmp2 * tmp4 / H0')

if __name__ == "__main__":
    q = EigenMatrix('ff=0.5.txt', 11, 11, np.array([1e-9, 0]))

    q.save_matrix_to_file()
