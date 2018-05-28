from src.eig_problem.EigenMatrix import EigenMatrix
from math import hypot
import numpy as np
from math import cos, sin
from scipy.linalg import eig
import sys
from src.eig_problem.EigenMatrix1D import EigenMatrix1D
from src.io.DataReader import ParsingData
import os.path

scriptpath = os.path.dirname(__file__)


class EigenValueProblem:
    def __init__(self, input_parameters, mat_1, mat_2):
        if isinstance(input_parameters, str):
            self.parameters = ParsingData(input_parameters)
        elif isinstance(input_parameters, dict):
            self.parameters = ParsingData(input_parameters)
        else:
            self.parameters = input_parameters
        self.mat_1 = mat_1
        self.mat_2 = mat_2

    def calculate_dispersion(self):
        pass

    def solve_eigen_problem(self, bloch_vector, param, bloch_vector_perp=0):
        pass

    def list_bloch_vector(self):
        pass

    def calculate_eigen_frequency(self, bloch_vector, bloch_vector_perp=0):
        gamma, mu0H0 = self.parameters.physical_constant()
        eigen_vector = self.solve_eigen_problem(bloch_vector, param=False, bloch_vector_perp=bloch_vector_perp)
        eigen_value = [i.imag * gamma * mu0H0 / 2.0 / np.pi for i in eigen_vector if i.imag > 0]
        return list(sorted(eigen_value)[:50])  # TODO: create smarter choice

    def calculate_eigen_vectors(self):
        eigen_value, eigen_vector = self.solve_eigen_problem(self.list_bloch_vector()[0], param=True)
        eigen_value_index = np.argsort(eigen_value.imag)
        eigen_vector = np.transpose(eigen_vector)
        eigen_vector = eigen_vector[eigen_value_index[len(eigen_value) // 2:]]
        return eigen_vector

    def print_eigen_vectors(self):
        np.savetxt(self.parameters.output_file(),
                   self.calculate_eigen_vectors().view(float),
                   header='Bloch wave vector, q=' + str(self.list_bloch_vector()[0]))


class EigenValueProblem2D(EigenValueProblem):
    def __init__(self, direction, input_parameters, mat_1, mat_2):
        EigenValueProblem.__init__(self, input_parameters, mat_1, mat_2)
        self.direction = direction
        if self.direction == 'x':
            self.coordinate = [0, 1]
        elif self.direction == 'y':
            self.coordinate = [1, 0]
        elif self.direction == 'xy':
            self.coordinate = [1, 1]
        else:
            sys.exit('Wrong argument for direction was set')

    # TODO: update symbol; add another paths

    def calculate_dispersion(self):
        data = []
        for k in self.list_bloch_vector():
            tmp = [hypot(k[0], k[1])]
            tmp.extend(self.calculate_eigen_frequency(k))
            data.append(tmp)
        return data

    def solve_eigen_problem(self, bloch_vector, param, bloch_vector_perp=0):
        eigen_matrix = EigenMatrix(EigenMatrix.ReciprocalVectorGrid(*self.parameters.rec_vector()),
                                   bloch_vector, self.parameters, self.mat_1, self.mat_2).generate_and_fill_matrix()

        return eig(eigen_matrix, right=param)

    def list_bloch_vector(self):
        points = np.linspace(*self.parameters.bloch_vector())
        return 2 * np.pi * np.stack((points, points), axis=-1) * self.coordinate / \
               self.parameters.lattice_const()


class EigenValueProblem1D(EigenValueProblem):
    def __init__(self, input_parameters, mat_1, mat_2):
        EigenValueProblem.__init__(self, input_parameters, mat_1, mat_2)

    def calculate_dispersion(self):
        data = []
        for k in self.list_bloch_vector():
            tmp = [k]
            tmp.extend(self.calculate_eigen_frequency(k))
            data.append(tmp)
        return data

    def oblique_dispersion(self):
        angle, max_value = self.parameters.perpendicular_bloch_vector()
        min, max, step = self.parameters.bloch_vector()
        max_value = np.pi/self.parameters.lattice_const()[0] / cos(np.radians(angle))
        bloch_vec = np.linspace(2.278401604183777526e+02, max_value, step)
        print(max_value)
        x_bloch_vec = bloch_vec * cos(np.radians(angle))
        y_bloch_vec = bloch_vec * sin(np.radians(angle))
        data = []
        for k1, k2 in zip(x_bloch_vec, y_bloch_vec):
            tmp = [np.sqrt(k1**2+k2**2)]
            tmp.extend(self.calculate_eigen_frequency(k1, k2))
            data.append(tmp)
        return data

    def solve_eigen_problem(self, bloch_vector, param, bloch_vector_perp=0):
        return eig(EigenMatrix1D(bloch_vector, self.parameters, self.mat_1,
                                 self.mat_2, bloch_vector_perp).generate_and_fill_matrix(),
                   right=param)  # trzeba pamiętać o włączeniu/wyłączeniu generowania wektorów

    def list_bloch_vector(self):
        return [2 * np.pi * k / self.parameters.lattice_const()[0]
                for k in np.linspace(*self.parameters.bloch_vector())]


if __name__ == "__main__":
    pass
