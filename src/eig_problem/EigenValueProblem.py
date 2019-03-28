from src.eig_problem.EigenMatrix import EigenMatrix
import numpy as np
from math import cos, sin
from scipy.linalg import eig
import sys
from src.eig_problem.EigenMatrix1D import EigenMatrix1D
from src.io.DataReader import ParsingData
import os.path
from itertools import repeat
from src.utils.cProfiler import do_cprofile


scriptpath = os.path.dirname(__file__)


class EigenValueProblem:
    def __init__(self, input_parameters, mat_1, mat_2):
        if isinstance(input_parameters, (str, dict)):
            self.parameters = ParsingData(input_parameters)
        elif isinstance(input_parameters, ParsingData):
            self.parameters = input_parameters
        else:
            raise IOError('Input parameters not understood: ' + str(input_parameters))
        self.mat_1 = mat_1
        self.mat_2 = mat_2

    def calculate_dispersion_along_direction(self):
        pass

    def solve_eigen_problem(self, bloch_vector, param, bloch_vector_perp=0):
        pass

    def list_bloch_vector(self):
        pass

    def calculate_eigen_frequency(self, bloch_vector, bloch_vector_perp=0):
        gamma, mu0H0 = self.parameters.physical_constant()
        eigen_vector = self.solve_eigen_problem(bloch_vector, param=False, bloch_vector_perp=bloch_vector_perp)
        eigen_value = [i.imag * gamma * mu0H0 / 2.0 / np.pi for i in eigen_vector if i.imag > 0]
        return np.array(list(sorted(eigen_value)[:50]))  # TODO: create smarter choice

    def calculate_eigen_vectors(self, bloch_vector=np.array([1, 1])):
        eigen_value, eigen_vector = self.solve_eigen_problem(bloch_vector, param=True)
        eigen_value_index = np.argsort(eigen_value.imag)
        eigen_vector = np.transpose(eigen_vector)
        eigen_vector = eigen_vector[eigen_value_index[len(eigen_value) // 2:]]
        return eigen_vector

    @do_cprofile
    def calculate_eigen_vectors_and_frequency(self, bloch_vector=np.array([1, 1])):
        gamma, mu0H0 = self.parameters.physical_constant()

        eigen_value, eigen_vector = self.solve_eigen_problem(bloch_vector, param=True)
        frequencies = [i.imag * gamma * mu0H0 / 2.0 / np.pi for i in eigen_value if i.imag > 0]
        eigen_value_index = np.argsort(eigen_value.imag)
        eigen_vector = np.transpose(eigen_vector)
        eigen_vector = eigen_vector[eigen_value_index[len(eigen_value) // 2:]]
        return eigen_vector[:40, :], np.array(list(sorted(frequencies)[:40]))

    def print_eigen_vectors(self):
        np.savetxt(self.parameters.output_file(),
                   self.calculate_eigen_vectors().view(float),
                   header='Bloch wave vector, q=' + str(self.list_bloch_vector()[0]))


class EigenValueProblem2D(EigenValueProblem):
    def __init__(self, direction, input_parameters, mat_1, mat_2):
        EigenValueProblem.__init__(self, input_parameters, mat_1, mat_2)
        self.direction = direction
        # TODO: this concept seems to be wrong
        if self.direction == 'x':
            self.coordinate = [0, 1]
        elif self.direction == 'y':
            self.coordinate = [1, 0]
        elif self.direction == 'xy':
            self.coordinate = [1, 1]
        elif self.direction == 'oblique':
            angle = self.parameters.perpendicular_bloch_vector()[0]
            self.coordinate = [cos(np.radians(angle)), sin(np.radians(angle))]

        else:
            sys.exit('Wrong argument for direction was set')

    def calculate_dispersion_map(self):
        points = np.linspace(*self.parameters.bloch_vector())
        return self.calculate_dispersion(2 * np.pi * cartesian_product_transpose_pp([points, points])/ \
               self.parameters.lattice_const())

    def calculate_dispersion_along_direction(self):
        return self.calculate_dispersion(self.list_bloch_vector())

    def calculate_dispersion(self, vectors):
        first_record = self.calculate_eigen_frequency(vectors[0])
        data = np.zeros((vectors.shape[0], len(first_record) + 2))
        data[0, :2] = vectors[0]
        data[0, 2:] = first_record
        for ind, k in enumerate(vectors[1:, :]):
            data[ind + 1, 0:2] = k
            data[ind + 1, 2:] = self.calculate_eigen_frequency(k)
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

    def calculate_dispersion_along_direction(self):
        data = []
        for k in self.list_bloch_vector():
            tmp = [k, 0]
            tmp.extend(self.calculate_eigen_frequency(k))
            data.append(tmp)
        return np.array(data)

    def oblique_dispersion(self):
        '''propagation along given direction'''
        angle, max_value = self.parameters.perpendicular_bloch_vector()
        min, max, step = self.parameters.bloch_vector()
        max_value = np.pi/self.parameters.lattice_const()[0] / cos(np.radians(angle))
        bloch_vec = np.linspace(100, 13*max_value, step)
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


def cartesian_product_transpose_pp(arrays):
    la = len(arrays)
    dtype = np.result_type(*arrays)
    arr = np.empty((la, *map(len, arrays)), dtype=dtype)
    idx = slice(None), *repeat(None, la)
    for i, a in enumerate(arrays):
        arr[i, ...] = a[idx[:la-i]]
    return arr.reshape(la, -1).T


if __name__ == "__main__":
    pass
