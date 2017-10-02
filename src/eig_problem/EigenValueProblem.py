
from src.eig_problem.EigenMatrix import EigenMatrix
from math import hypot
import numpy as np
from scipy.linalg import eig
import sys
from src.eig_problem.EigenMatrix1D import EigenMatrix1D
from src.eig_problem.InputParameter import InputParameter
import os.path

scriptpath = os.path.dirname(__file__)

class EigenValueProblem:
    def __init__(self, number_of_dispersion_point, a, gamma, mu0H0, input_fft_file):
        self.number_of_dispersion_point = number_of_dispersion_point
        self.gamma = gamma
        self.mu0H0 = mu0H0

        self.input_fft_file = input_fft_file
        self.a = a
        self.start_vec_q = 0.01
        self.end_vec_q = 0.5

    def eigen_frequency_for_vectors_q(self):
        pass

    def solve_eigen_problem(self, wektor_q, param):
        pass

    def list_vector_q(self):
        pass

    def calculate_eigen_frequency(self, wektor_q):
        eigen_vector = self.solve_eigen_problem(wektor_q, param=False)
        eigen_value = [i.imag * self.gamma * self.mu0H0 / 2.0 / np.pi for i in eigen_vector if i.imag > 0]
        return list(sorted(eigen_value)[:50]) # TODO: create smarter choice

    def calculate_eigen_vectors(self):
        eigen_value, eigen_vector = self.solve_eigen_problem(self.list_vector_q()[0], param=True)
        eigen_value_index = np.argsort(eigen_value.imag)
        eigen_vector = np.transpose(eigen_vector)
        eigen_vector = eigen_vector[eigen_value_index[len(eigen_value) // 2:]]
        return eigen_vector.view(float)

    def print_eigen_vectors(self):
        np.savetxt(str(self.list_vector_q()[0]) + '.', eigen_vector.view(float))


class EigenValueProblem2D(EigenValueProblem):
    def __init__(self, number_of_dispersion_point, direction, a=InputParameter.a, b=InputParameter.b,
                 gamma=InputParameter.gamma, mu0H0=InputParameter.mu0H0, input_fft_file=InputParameter.fft_file):

        EigenValueProblem.__init__(self, number_of_dispersion_point, a, gamma, mu0H0, input_fft_file)

        self.b = b
        self.direction = direction
        self.input_fft_file = 'ff=0.5.txt'
        if self.direction == 'x':
            self.coordinate = [0, 1]
        elif self.direction == 'y':
            self.coordinate = [1, 0]
        elif self.direction == 'xy':
            self.coordinate = [1, 1]
        else:
            sys.exit('Wrong argument for direction was set')
    # TODO: update symbol; add another paths

    def eigen_frequency_for_vectors_q(self):
        data = []
        for k in self.list_vector_q():
            tmp = [hypot(k[0], k[1])]
            tmp.extend(self.calculate_eigen_frequency(k))
            data.append(tmp)
        return data

    def solve_eigen_problem(self, wektor_q, param):
        eigen_matrix = EigenMatrix(self.input_fft_file, EigenMatrix.ReciprocalVectorGrid(9, 9),
                                   wektor_q).generate_and_fill_matrix()
        return eig(eigen_matrix, right=param)

    def list_vector_q(self):
        points = np.linspace(self.start_vec_q, self.end_vec_q, self.number_of_dispersion_point)
        return 2 * np.pi * np.stack((points, points), axis=-1) * self.coordinate / [self.a, self.b]


class EigenValueProblem1D(EigenValueProblem):
    def __init__(self, number_of_dispersion_point, input_fft_file, output_file, mu0H0=InputParameter.mu0H0,
                  a=InputParameter.a, gamma=InputParameter.gamma, angle=InputParameter.angle):

        EigenValueProblem.__init__(self, number_of_dispersion_point, a, gamma, mu0H0, input_fft_file)

        self.input_fft_file = os.path.join(scriptpath, input_fft_file)
        self.output_file = output_file
        self.angle = angle

    def eigen_frequency_for_vectors_q(self):
        data = []
        for k in self.list_vector_q():
            tmp = [k]
            tmp.extend(self.calculate_eigen_frequency(k))
            data.append(tmp)
        return data

    def solve_eigen_problem(self, wektor_q, param):
        macierz_m = EigenMatrix1D(self.input_fft_file, wektor_q,
                                  angle=self.angle).matrix_angle_dependence(wektor_q)
        return eig(macierz_m, right=param)  # trzeba pamiętać o włączeniu/wyłączeniu generowania wektorów

    def list_vector_q(self):
        return [2 * np.pi * k / self.a for k in np.linspace(0.0001, 0.9999, self.number_of_dispersion_point)]


if __name__ == "__main__":
    EigenValueProblem1D(1, 'c_coef_100.txt', 'dys.txt').calculate_eigen_vectors()