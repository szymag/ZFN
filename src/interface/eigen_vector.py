import numpy as np
from src.fft_from_image.FFT import FFT
from src.eig_problem.EigenValueProblem import EigenValueProblem2D
from src.eig_problem.InputParameter import InputParameter
from src.drawing.Plot import Plot

file_name = 'ff=0.5.txt'
number_of_dispersion_point = 20
direction = 'xy'
number_of_dispersion_branch = 5
x_lim = None
y_lim = None
show_plot = True

def do_program():
    if file_name[-3:] == 'png':
        fft = FFT().zwroc_tablice(file_name)
    elif file_name[-3:] == 'txt':
        fft = np.loadtxt(file_name)
    elif file_name[-3:] == 'dat':
        return Plot(number_of_dispersion_branch, x_lim, y_lim,
                    name_of_output_file=InputParameter.output_file).dispersion_relation(file_name)
    else:
        raise ValueError

    eig_vec = EigenValueProblem2D(1, direction,
                                  input_fft_file=fft).calculate_eigen_vectors()
    np.savetxt('tst.txt', eig_vec)
    if show_plot:
        return 0
    else:
        return 0

do_program()
