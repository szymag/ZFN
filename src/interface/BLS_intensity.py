import numpy as np
from src.eig_problem.EigenValueProblem import EigenValueProblem2D, cartesian_product_transpose_pp
from src.io.DataReader import ParsingData
from src.drawing.Plot import Plot


input_parameters = ParsingData('./src/interface/set5.yaml')
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
        eig_vec, eig_freq = tmp.calculate_eigen_vectors_and_frequency(i)
        dispersion[ind, :2] = i
        weights[ind, :2] = i

        dispersion[ind, 2:] = eig_freq
        weights[ind, 2:] = frequency_weight(eig_vec)[:40]
    return dispersion, weights


def BLS_intensity_map():
    tmp = EigenValueProblem2D(direction, input_parameters, 'Py', 'Co')
    points = np.linspace(*tmp.parameters.bloch_vector())
    bloch_vectors = 2 * np.pi * cartesian_product_transpose_pp([points, points])/ \
               tmp.parameters.lattice_const()

    dispersion = np.zeros((bloch_vectors.shape[0], 7))
    weights = np.zeros((bloch_vectors.shape[0], 7))
    for ind, i in enumerate(bloch_vectors):
        eig_vec, eig_freq = tmp.calculate_eigen_vectors_and_frequency(i)

        dispersion[ind, :2] = i
        weights[ind, :2] = i

        bls_weight = frequency_weight(eig_vec)
        weights_ind = np.argsort(bls_weight)[-1:-6:-1]
        dispersion[ind, 2:] = eig_freq[weights_ind]
        weights[ind, 2:] = bls_weight[weights_ind]
    return dispersion, weights


def visualize(frequencies, weights, angle=None):
    return Plot(number_of_dispersion_branch,
                y_lim=[5, 25],
                x_lim=[0, 3.5]).bls(frequencies, weights, angle)


def visualize_map(frequencies, weights, frequency):
    return Plot(number_of_dispersion_branch, x_lim=[-1.05, 1.05], y_lim=[-1.05, 1.05]).bls_for_given_frequency(
        frequencies, weights, frequency)


if __name__ == "__main__":
    pass
