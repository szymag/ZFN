from src.eig_problem.EigenMatrix import EigenMatrix
from math import hypot
import numpy as np
from scipy.linalg import eig
import sys
from src.eig_problem.EigenMatrix1D import EigenMatrix1D
from src.eig_problem.InputParameter import InputParameter
from src.io.DataReader import load_yaml_file
import os.path

scriptpath = os.path.dirname(__file__)


class EigenValueProblem:
    def __init__(self, input_parameters):
        self.parameters = load_yaml_file(input_parameters)
        self.number_of_dispersion_point = self.parameters['q_vector']['dispersion_count']
        self.gamma = self.parameters['physical_parameters']['gamma']
        self.mu0H0 = self.parameters['physical_parameters']['mu0H0']

        self.input_fft_file = self.parameters['numerical_parameters']['fft_file']
        if self.parameters['numerical_parameters']['output_file'][-4:] == '.vec':
            self.output_file_name = self.parameters['numerical_parameters']['output_file']
        else:
            self.output_file_name = self.parameters['numerical_parameters']['output_file'] + '.vec'
        self.rec_vector_x = self.parameters['numerical_parameters']['rec_vector_x']
        self.rec_vector_y = self.parameters['numerical_parameters']['rec_vector_y']
        self.lattice_const_x = self.parameters['system_dimensions']['a']
        self.lattice_const_y = self.parameters['system_dimensions']['b']
        self.start_vec_q = self.parameters['q_vector']['start']
        self.end_vec_q = self.parameters['q_vector']['end']

    def eigen_frequency_for_bloch_vectors(self):
        pass

    def solve_eigen_problem(self, wektor_q, param):
        pass

    def list_vector_q(self):
        pass

    def calculate_eigen_frequency(self, bloch_vector):
        eigen_vector = self.solve_eigen_problem(bloch_vector, param=False)
        eigen_value = [i.imag * self.gamma * self.mu0H0 / 2.0 / np.pi for i in eigen_vector if i.imag > 0]
        return list(sorted(eigen_value)[:50])  # TODO: create smarter choice

    def calculate_eigen_vectors(self):
        eigen_value, eigen_vector = self.solve_eigen_problem(self.list_vector_q()[0], param=True)
        eigen_value_index = np.argsort(eigen_value.imag)
        eigen_vector = np.transpose(eigen_vector)
        eigen_vector = eigen_vector[eigen_value_index[len(eigen_value) // 2:]]
        return eigen_vector

    def print_eigen_vectors(self):
        np.savetxt(self.output_file_name,
                   self.calculate_eigen_vectors().view(float),
                   header='Bloch wave vector, q=' + str(self.list_vector_q()[0]))


class EigenValueProblem2D(EigenValueProblem):
    def __init__(self, direction, input_parameters):
        EigenValueProblem.__init__(self, input_parameters)
        self.direction = direction
        self.input_fft_file = 'ff=0.5.fft'
        if self.direction == 'x':
            self.coordinate = [0, 1]
        elif self.direction == 'y':
            self.coordinate = [1, 0]
        elif self.direction == 'xy':
            self.coordinate = [1, 1]
        else:
            sys.exit('Wrong argument for direction was set')

    # TODO: update symbol; add another paths

    def eigen_frequency_for_bloch_vectors(self):
        data = []
        for k in self.list_vector_q():
            tmp = [hypot(k[0], k[1])]
            tmp.extend(self.calculate_eigen_frequency(k))
            data.append(tmp)
        return data

    def solve_eigen_problem(self, wektor_q, param):
        print(self.parameters)
        eigen_matrix = EigenMatrix(EigenMatrix.ReciprocalVectorGrid(
            self.rec_vector_x,
            self.rec_vector_y), wektor_q, self.parameters, 'Fe', 'Ni').generate_and_fill_matrix()
        return eig(eigen_matrix, right=param)

    def list_vector_q(self):
        points = np.linspace(self.start_vec_q, self.end_vec_q, self.number_of_dispersion_point)
        return 2 * np.pi * np.stack((points, points), axis=-1) * self.coordinate / \
               [self.lattice_const_x, self.lattice_const_y]


class EigenValueProblem1D(EigenValueProblem):
    def __init__(self, number_of_dispersion_point, input_fft_file, output_file, input_parameters,
                 angle=InputParameter.angle):
        EigenValueProblem.__init__(self, number_of_dispersion_point, input_parameters)

        self.input_fft_file = os.path.join(scriptpath, input_fft_file)
        self.output_file = output_file
        self.angle = angle

    def eigen_frequency_for_bloch_vectors(self):
        data = []
        for k in self.list_vector_q():
            tmp = [k]
            tmp.extend(self.calculate_eigen_frequency(k))
            data.append(tmp)
        return data

    def solve_eigen_problem(self, bloch_vector, param):
        return eig(EigenMatrix1D(self.input_fft_file, bloch_vector,
                                 angle=self.angle).matrix_angle_dependence(bloch_vector)
                   , right=param)  # trzeba pamiętać o włączeniu/wyłączeniu generowania wektorów

    def list_vector_q(self):
        return [2 * np.pi * k / self.lattice_const_x
                for k in np.linspace(0.0001, 0.9999, self.number_of_dispersion_point)]


if __name__ == "__main__":
    EigenValueProblem2D('x', 'InputParameter.yaml').print_eigen_vectors()
