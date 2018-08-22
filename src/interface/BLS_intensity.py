import numpy as np
from src.eig_problem.EigenValueProblem import EigenValueProblem2D, EigenValueProblem1D
from src.io.DataReader import ParsingData
from src.drawing.Plot import Plot


input_parameters = ParsingData('JEMS.yaml')
direction = 'oblique'
number_of_dispersion_branch = 40
x_lim = None
y_lim = None
show_plot = True


def frequency_weight(vectors):
    return abs(vectors[:, 3*vectors.shape[1]//4])**2


def BLS_intensity():
    # TODO: Define class responsible for handle formats
    tmp = EigenValueProblem2D(direction, input_parameters, 'Py', 'Co')
    bloch_vectors = tmp.list_bloch_vector()
    dispersion = np.zeros((input_parameters.bloch_vector()[2], 42))
    weights = np.zeros((input_parameters.bloch_vector()[2], 42))
    for ind, i in enumerate(bloch_vectors):
        eig = tmp.calculate_eigen_vectors_and_frequency(i)
        dispersion[ind, :2] = i
        weights[ind, :2] = i

        dispersion[ind, 2:] = eig[1]
        weights[ind, 2:] = frequency_weight(eig[0])[:40]

    return dispersion, weights


def visualize():
    return Plot(number_of_dispersion_branch, y_lim=[0, 25e9]).bls(*BLS_intensity())


if __name__ == "__main__":
    visualize()