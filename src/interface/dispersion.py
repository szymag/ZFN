import numpy as np
from src.fft_from_image.FFT import FFT
from src.eig_problem.EigenValueProblem import EigenValueProblem2D, EigenValueProblem1D
from src.io.DataReader import ParsingData
from src.drawing.Plot import Plot

input_parameters = ParsingData('alongwires.yaml')
file_name = 'tst.dat'
direction = 'xy'
number_of_dispersion_branch = 3
x_lim = None
y_lim = None
show_plot = False


def do_program():  # TODO: file type should represent containing data
    if file_name[-3:] == 'png':
        fft = FFT().zwroc_tablice(file_name)
    elif file_name[-3:] == 'txt':
        fft = np.loadtxt(file_name)
    elif file_name[-3:] == 'dat':
        return Plot(number_of_dispersion_branch, x_lim, y_lim,
                    name_of_output_file=input_parameters.output_file()).dispersion_relation(file_name)
    else:
        raise ValueError

    eig_freq = EigenValueProblem2D(direction, input_parameters, 'Co', 'Py').calculate_dispersion()
    np.savetxt('tst.txt', eig_freq)
    if show_plot:
        return Plot(number_of_dispersion_branch, x_lim, y_lim,
                    name_of_output_file=input_parameters.output_file()).dispersion_relation(eig_freq)
    else:
        return 0


def do_program_1D():
    eig_freq = EigenValueProblem1D(input_parameters, 'Py', 'Py').calculate_dispersion()
    np.savetxt(file_name, eig_freq)
    if show_plot:
        return Plot(number_of_dispersion_branch, x_lim, y_lim,
                    name_of_output_file=input_parameters.output_file()).dispersion_relation(eig_freq)
    else:
        return 0


if __name__ == "__main__":
    eig_freq = EigenValueProblem1D(input_parameters, 'Py', 'Py').oblique_dispersion()
    #print(eig_freq)
    np.savetxt(file_name, eig_freq)
    Plot(2).dispersion_relation('tst.dat')
    #do_program_1D()