# -*- coding: utf-8 -*-
import numpy as np

from src.eig_problem.LoadFFT import LoadFFT1D
from src.eig_problem.ReciprocalVector import ReciprocalVector
from src.eig_problem.InputParameter import InputParameter

class EigenMatrix1D:

    def __init__(self, input_fft, wektor_q, H0=InputParameter.H0, a=InputParameter.a,
                 MoA=InputParameter.MoA,
                 MoB=InputParameter.MoB, lA=InputParameter.lA,
                 lB=InputParameter.lB, d=InputParameter.d,
                x=InputParameter.x, angle=InputParameter.angle):

        self.a = a
        self.d = d
        self.x = x
        self.H0 = H0
        self.tmp = LoadFFT1D(input_fft)
        self.ilosc_wektorow = self.tmp.ilosc_wektorow
        self.macierz_M = np.zeros((2 * self.ilosc_wektorow, 2 * self.ilosc_wektorow), dtype=complex)
        self.magnetyzacja = self.tmp.fourier_coefficient(MoA, MoB)
        self.dlugosc_wymiany = self.tmp.fourier_coefficient(lA, lB)
        self.lista_wektorow = ReciprocalVector(self.ilosc_wektorow).lista_wektorow1d('min')
        self.shift = len(self.lista_wektorow) - 1
        self.angle = angle

    def funkcja_c(self, wektor_1, wektor_2):
        return np.exp(-abs((wektor_1 + wektor_2)) * self.d / 2)

    def delta_kroneckera(self):
        for i in range(self.ilosc_wektorow, 2 * self.ilosc_wektorow):
            self.macierz_M[i - self.ilosc_wektorow][i] += 1.
            self.macierz_M[i][i - self.ilosc_wektorow] -= 1.

    def exchange_field(self, wektor_1, wektor_2, wektor_q):
        vec_l = np.transpose(np.broadcast_to(self.lista_wektorow ,
         (self.ilosc_wektorow, self.ilosc_wektorow)))
        tmp1 = self.dlugosc_wymiany[vec_l - wektor_2 + self.shift]
        tmp2 = self.magnetyzacja[wektor_1 - vec_l + self.shift]
        tmp3 = (wektor_q + 2 * np.pi * wektor_2 / self.a) * (2 * np.pi * vec_l / self.a + wektor_q)
        return np.sum(tmp1 * tmp2 * tmp3, axis=0) / self.H0

    def dynamic_demagnetizing_field_in_plane(self, wektor_1, wektor_2, wektor_q):
        tmp3 = self.magnetyzacja[wektor_1 - wektor_2 + self.shift]
        tmp2 = 1 - self.funkcja_c(wektor_q, (2 * np.pi * wektor_2 / self.a))
        return tmp3 * tmp2 / self.H0

    def dynamic_demagnetizing_field_out_of_plane(self, wektor_1, wektor_2, wektor_q):
        tmp3 = self.magnetyzacja[wektor_1 - wektor_2 + self.shift]
        tmp2 = self.funkcja_c(wektor_q, (2 * np.pi * wektor_2 / self.a))
        return tmp2 * tmp3 / self.H0

    def static_demagnetizing_field(self, wektor_1, wektor_2):

        tmp1 = self.magnetyzacja[wektor_1 - wektor_2 + self.shift]
        tmp2 = 1 - np.exp(-abs((2 * np.pi * wektor_1 / self.a - 2 * np.pi * wektor_2 / self.a)) * self.d / 2)
        return tmp2 * tmp1 / self.H0

    def matrix_angle_dependence(self, wektor_q):
        indeks = self.ilosc_wektorow
        self.delta_kroneckera()
        for i in range(indeks, 2 * indeks):
            w1 = self.lista_wektorow[i - indeks]
            w2 = self.lista_wektorow
            ex = self.exchange_field(w1, w2, wektor_q)
            dyn_in_plane = self.dynamic_demagnetizing_field_in_plane(w1, w2, wektor_q)
            dyn_out_plane = self.dynamic_demagnetizing_field_out_of_plane(w1, w2, wektor_q)
            static = self.static_demagnetizing_field(w1, w2)
            self.macierz_M[i][np.arange(indeks)] += -ex - dyn_in_plane + static * np.cos\
                (np.pi * self.angle / 180)**2  # yx
            self.macierz_M[i - indeks][np.arange(indeks, 2 * indeks)] += ex + dyn_out_plane * np.sin\
                (np.pi * self.angle / 180)**2 - static * np.cos(np.pi * self.angle / 180)**2  # xy
        return self.macierz_M

    def wypisz_macierz(self):

        self.matrix_angle_dependence(1e-9)
        np.savetxt('macierz.txt', np.array(self.macierz_M))


if __name__ == "__main__":
    q = EigenMatrix1D('p_coef_100*2.txt', 1e-9)
    q.wypisz_macierz()
