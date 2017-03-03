from src.eig_problem.cProfiler import do_cprofile
import numpy as np
from math import sqrt
from scipy.linalg import eig
from src.eig_problem.EigenMatrix import EigenMatrix
from src.eig_problem.InputParameter import InputParameter
import sys


class EigenValueProblem:

    def __init__(self, ilosc_wektorow_q, direction, a=InputParameter.a, b=InputParameter.b,
                 gamma=InputParameter.gamma, mu0H0=InputParameter.mu0H0):
        self.number_of_vec_q = ilosc_wektorow_q
        self.gamma = gamma
        self.mu0H0 = mu0H0
        self.direction = direction
        self.a = a
        self.b = b
        self.start_vec_q = 0.01
        self.end_vec_q = 0.5

    def save_frequency_to_file(self):
        plik = []
        if self.direction == 'x':
            for k in self.list_vector_q():
                tmp = [k[1]]
                tmp.extend(self.calculate_eigen_frequency(k))
                plik.append(tmp)
        elif self.direction == 'y':
            for k in self.list_vector_q():
                tmp = [k[0]]
                tmp.extend(self.calculate_eigen_frequency(k))
                plik.append(tmp)
        elif self.direction == 'xy':
            for k in self.list_vector_q():
                tmp = [sqrt(k[0] ** 2 + k[1] ** 2)]
                tmp.extend(self.calculate_eigen_frequency(k))
                plik.append(tmp)
        np.savetxt('test.txt', plik)

    def calculate_eigen_frequency(self, wektor_q):
        assert len(wektor_q) == 2, \
            'form of wektor_q is forbidden. wektor_q should have two arguments'
        wartosci_wlasne = self.solve_eigen_problem(wektor_q, param=False)
        czestosci_wlasne = [i.imag * self.gamma * self.mu0H0 / 2.0 / np.pi for i in wartosci_wlasne if i.imag > 0]
        return list(sorted(czestosci_wlasne)[:50])

    def calculate_eigen_vectors(self):
        assert len(self.list_vector_q()) == 1, 'Eigenvector should be calculated for only one position vector'
        wartosci_wlasne, wektory_wlasne = self.solve_eigen_problem(self.list_vector_q()[0], param=True)
        wartosci_wlasne_index = np.argsort(wartosci_wlasne.imag)
        wektory_wlasne = np.transpose(wektory_wlasne)
        wektory_wlasne = wektory_wlasne[wartosci_wlasne_index[self.ilosc_wektorow:]]
        return np.savetxt(str(self.list_vector_q()[0]) + '.', wektory_wlasne.view(float))

    @do_cprofile
    def solve_eigen_problem(self, wektor_q, param):
        macierz_m = EigenMatrix("ff=0.5.txt", 11, 11, wektor_q).generate_and_fill_matrix()
        return eig(macierz_m, right=param)  # trzeba pamiętać o włączeniu/wyłączeniu generowania wektorów

    def list_vector_q(self):
        if self.direction == 'y':
            return [[(2 * np.pi * k / self.a), 0]
                    for k in np.linspace(self.start_vec_q, self.end_vec_q, self.number_of_vec_q)]
        elif self.direction == 'x':
            return [[0, (2 * np.pi * k / self.b)]
                    for k in np.linspace(self.start_vec_q, self.end_vec_q, self.number_of_vec_q)]
        elif self.direction == 'xy':
            return [[(2 * np.pi * k / self.a), (2 * np.pi * k / self.b)]
                    for k in np.linspace(self.start_vec_q, self.end_vec_q, self.number_of_vec_q)]
        else:
            sys.exit('Wrong argumnet for direction was set')


if __name__ == "__main__":
    def start():
        return EigenValueProblem(20, 'x').save_frequency_to_file()
    
    start()
