# -*- coding: utf-8 -*-
import numpy as np
from math import radians

from src.eig_problem.LoadFFT import LoadFFT1D
from src.eig_problem.ReciprocalVector import ReciprocalVector
from src.io.DataReader import ParsingData


class EigenMatrix1D:
    def __init__(self, input_parameters, material_A, material_B):
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
        self.ilosc_wektorow = self.tmp.ilosc_wektorow
        self.macierz_M = np.zeros((2 * self.ilosc_wektorow, 2 * self.ilosc_wektorow), dtype=complex)
        self.magnetization_sat = self.tmp.fourier_coefficient(self.material_A['Mo'], self.material_B['Mo'])
        self.exchange_len = self.tmp.fourier_coefficient(self.material_A['l'], self.material_B['l'])
        self.reciprocal_vec = ReciprocalVector(self.ilosc_wektorow).lista_wektorow1d('min')
        self.shift = len(self.reciprocal_vec) - 1

    def funkcja_c(self, wektor_1, wektor_2):
        return np.exp(-abs((wektor_1 + wektor_2)) * self.parameters.thickness() / 2)

    def delta_kroneckera(self):
        for i in range(self.ilosc_wektorow, 2 * self.ilosc_wektorow):
            self.macierz_M[i - self.ilosc_wektorow][i] += 1.
            self.macierz_M[i][i - self.ilosc_wektorow] -= 1.

    def exchange_field(self, wektor_1, wektor_2, wektor_q):
        vec_l = np.transpose(np.broadcast_to(self.reciprocal_vec, (self.ilosc_wektorow, self.ilosc_wektorow)))
        tmp1 = self.exchange_len[vec_l - wektor_2 + self.shift]
        tmp2 = self.magnetization_sat[wektor_1 - vec_l + self.shift]
        tmp3 = (wektor_q + 2 * np.pi * wektor_2 /
                self.parameters.lattice_const()[0]) * (2 * np.pi * vec_l /
                                                       self.parameters.lattice_const()[0] + wektor_q)
        return np.sum(tmp1 * tmp2 * tmp3, axis=0) / self.H0

    def dynamic_demagnetizing_field_in_plane(self, wektor_1, wektor_2, wektor_q):
        tmp3 = self.magnetization_sat[wektor_1 - wektor_2 + self.shift]
        tmp2 = 1 - self.funkcja_c(wektor_q, (2 * np.pi * wektor_2 / self.parameters.lattice_const()[0]))
        return tmp3 * tmp2 / self.H0

    def dynamic_demagnetizing_field_out_of_plane(self, wektor_1, wektor_2, wektor_q):
        tmp3 = self.magnetization_sat[wektor_1 - wektor_2 + self.shift]
        tmp2 = self.funkcja_c(wektor_q, (2 * np.pi * wektor_2 / self.parameters.lattice_const()[0]))
        return tmp2 * tmp3 / self.H0

    def static_demagnetizing_field(self, wektor_1, wektor_2):
        tmp1 = self.magnetization_sat[wektor_1 - wektor_2 + self.shift]
        co = np.cosh(abs((2 * np.pi * wektor_1 / self.parameters.lattice_const()[0] - 2 * np.pi * wektor_2 /
                          self.parameters.lattice_const()[0])) * self.parameters.x())
        tmp2 = 1 - co*np.exp(-abs((2 * np.pi * wektor_1 / self.parameters.lattice_const()[0] - 2 * np.pi * wektor_2 /
                                   self.parameters.lattice_const()[0])) * self.parameters.thickness() / 2)
        tmp1[self.shift] = 0
        return tmp2 * tmp1 / self.H0

    def matrix_angle_dependence(self, wektor_q):
        # TODO: Inconsistency in storing Bloch vector
        indeks = self.ilosc_wektorow
        self.delta_kroneckera()
        for i in range(indeks, 2 * indeks):
            w1 = self.reciprocal_vec[i - indeks]
            w2 = self.reciprocal_vec
            ex = self.exchange_field(w1, w2, wektor_q)
            dyn_in_plane = self.dynamic_demagnetizing_field_in_plane(w1, w2, wektor_q)
            dyn_out_plane = self.dynamic_demagnetizing_field_out_of_plane(w1, w2, wektor_q)
            static = self.static_demagnetizing_field(w1, w2)
            self.macierz_M[i][np.arange(indeks)] += -ex - dyn_in_plane * np.sin\
                (radians(self.parameters.angle()))**2 + static * np.cos(radians(self.parameters.angle()))**2  # yx
            self.macierz_M[i - indeks][np.arange(indeks, 2 * indeks)] += \
                ex + dyn_out_plane - static * np.cos(radians(self.parameters.angle()))**2  # xy
        return self.macierz_M

    def wypisz_macierz(self):
        self.matrix_angle_dependence(1e-9)
        np.savetxt('macierz.txt', np.array(self.macierz_M))


if __name__ == "__main__":
    q = EigenMatrix1D('soko.fft', 1e-9)
    q.wypisz_macierz()
