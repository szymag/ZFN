import numpy as np
from src.fft_from_image.FFT import FFT
from src.eig_problem.EigenValueProblem import EigenValueProblem2D, EigenValueProblem1D
from src.io.DataReader import ParsingData
from src.drawing.Plot import Plot
from matplotlib import pyplot as plt

input_parameters = ParsingData('Centala.yaml')
file_name = 'dys_2.dat'
direction = 'xy'
number_of_dispersion_branch = 5
x_lim = None
y_lim = None
show_plot = True


def do_program():  # TODO: file type should represent containing data
    if file_name[-3:] == 'png':
        fft = FFT().zwroc_tablice(file_name)
        np.savetxt('tmp.fft', fft)
        input_parameters.set_new_value('tmp.vec', 'numerical_parameters', 'fft_file')
    elif file_name[-3:] == 'fft':
        fft = np.loadtxt(file_name)
    elif file_name[-3:] == 'dys':
        return Plot(number_of_dispersion_branch, x_lim, y_lim).dispersion_relation(file_name)
    else:
        raise ValueError

    eig_freq = EigenValueProblem2D(direction, input_parameters, 'Co', 'Py').calculate_dispersion()
    np.savetxt(input_parameters.output_file() + '.dys', eig_freq)
    if show_plot:
        return Plot(number_of_dispersion_branch, x_lim, y_lim).dispersion_relation(eig_freq)
    else:
        return 0


def do_program_1D():
    eig_freq = EigenValueProblem1D(input_parameters, 'CoFeB_1', 'CoFeB_2').calculate_dispersion()
    np.savetxt(file_name, eig_freq)
    if show_plot:
        return Plot(number_of_dispersion_branch, x_lim, y_lim).dispersion_relation(file_name)
    else:
        return 0


def do_program_oblique():
    eig_freq = EigenValueProblem1D('alongwires.yaml', 'Co', 'Py').oblique_dispersion()
    np.savetxt(file_name, eig_freq)
    Plot(5, x_lim, y_lim).dispersion_relation(file_name)


if __name__ == "__main__":
    do_program_oblique()
