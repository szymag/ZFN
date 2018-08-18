import numpy as np
from src.fft_from_image.FFT import FFT
from src.eig_problem.EigenValueProblem import EigenValueProblem2D
from src.io.DataReader import ParsingData

from src.eig_problem.InputParameter import InputParameter
from src.drawing.Plot import Plot

input_parameters = ParsingData('Parameter_for_TheImpact.yaml')
number_of_dispersion_point = 20
direction = 'xy'
number_of_dispersion_branch = 5
x_lim = None
y_lim = None
show_plot = True
# TODO: input_fft_file is wrong name. This argument is broader


def do_program():
    if input_parameters.input_fft_file()[-3:] == 'png':
        fft = FFT().zwroc_tablice(input_parameters.input_fft_file())
        np.savetxt('tmp.fft', fft)
        input_parameters.set_new_value('tmp.vec', 'numerical_parameters', 'fft_file')
    elif input_parameters.input_fft_file()[-3:] == 'fft':
        pass
    elif input_parameters.input_fft_file()[-3:] == 'vec':
        pass
        # TODO: Should be printed out
    else:
        raise ValueError

    eig_vec = EigenValueProblem2D(direction, input_parameters, 'Fe', 'Ni').calculate_eigen_vectors().view(float)
    np.savetxt(input_parameters.output_file('vectors'), eig_vec)
    if show_plot:
        return 0
    else:
        return 0


if __name__ == "__main__":
    do_program()
