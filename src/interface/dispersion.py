import numpy as np
from src.fft_from_image.FFT import FFT
from src.eig_problem.EigenValueProblem import EigenValueProblem2D
from src.eig_problem.InputParameter import InputParameter
from src.drawing.Plot import Plot

file_name = 'tst.dat'
number_of_dispersion_point = 20
direction = 'xy'
number_of_dispersion_branch = 5
x_lim = None
y_lim = [6, 10]
show_plot = True


def do_program():  # TODO: file type should represent containing data
    if file_name[-3:] == 'png':
        fft = FFT().zwroc_tablice(file_name)
    elif file_name[-3:] == 'txt':
        fft = np.loadtxt(file_name)
    elif file_name[-3:] == 'dat':
        return Plot(number_of_dispersion_branch, x_lim, y_lim,
                    name_of_output_file=InputParameter.output_file).dispersion_relation(file_name)
    else:
        raise ValueError

    eig_freq = EigenValueProblem2D(number_of_dispersion_point, direction,
                                   input_fft_file=fft).eigen_frequency_for_vectors_q()
    np.savetxt('tst.txt', eig_freq)
    if show_plot:
        return Plot(number_of_dispersion_branch, x_lim, y_lim,
                    name_of_output_file=InputParameter.output_file).dispersion_relation(eig_freq)
    else:
        return 0

do_program()
