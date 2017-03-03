from src.eig_problem.cProfiler import do_cprofile
import numpy as np
from scipy.linalg import eig
from src.eig_problem.EigenMatrix import EigenMatrix
from src.eig_problem.InputParameter import InputParameter
import sys
from math import hypot

class EigenValueProblem:
    def __init__(self, number_of_dispersion_point, direction, a=InputParameter.a, b=InputParameter.b,
                 gamma=InputParameter.gamma, mu0H0=InputParameter.mu0H0):
        self.number_of_dispersion_point = number_of_dispersion_point
        self.gamma = gamma
        self.mu0H0 = mu0H0
        self.direction = direction
        self.a = a
        self.b = b
        self.start_vec_q = 0.01
        self.end_vec_q = 0.5

        if self.direction == 'x':
            self.coordinate = [0, 1]
        elif self.direction == 'y':
            self.coordinate = [1, 0]
        elif self.direction == 'xy':
            self.coordinate = [1, 1]
        else:
            sys.exit('Wrong argumnet for direction was set')

    def save_frequency_to_file(self):
        file = []
        for k in self.list_vector_q():
            tmp = [hypot(k[0], k[1])]
            tmp.extend(self.calculate_eigen_frequency(k))
            file.append(tmp)
        np.savetxt('test.txt', file)

    def calculate_eigen_frequency(self, wektor_q):
        assert len(wektor_q) == 2, \
            'form of wektor_q is forbidden. wektor_q should have two arguments'
        wartosci_wlasne = self.solve_eigen_problem(wektor_q, param=False)
        czestosci_wlasne = [i.imag * self.gamma * self.mu0H0 / 2.0 / np.pi for i in wartosci_wlasne if i.imag > 0]
        return list(sorted(czestosci_wlasne)[:50])

    def calculate_eigen_vectors(self):
        assert len(self.list_vector_q()) == 1, 'Eigenvector should be calculated for only one position vector'
        eigen_value, eigen_vector = self.solve_eigen_problem(self.list_vector_q()[0], param=True)
        eigen_value_index = np.argsort(eigen_value.imag)
        eigen_vector = np.transpose(eigen_vector)
        eigen_vector = eigen_vector[eigen_value_index[len(eigen_value) // 2:]]
        return np.savetxt(str(self.list_vector_q()[0]) + '.', eigen_vector.view(float))

    @do_cprofile
    def solve_eigen_problem(self, wektor_q, param):
        macierz_m = EigenMatrix("ff=0.5.txt", 9, 9, wektor_q).generate_and_fill_matrix()
        return eig(macierz_m, right=param)

    def list_vector_q(self):
        points = np.linspace(self.start_vec_q, self.end_vec_q, self.number_of_dispersion_point)
        return 2 * np.pi * np.stack((points, points), axis=-1) * self.coordinate / [self.a, self.b]


if __name__ == "__main__":
    def start():
        return EigenValueProblem(20, 'y').save_frequency_to_file()
    
    start()
