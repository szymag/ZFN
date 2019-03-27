import numpy as np
import matplotlib.pyplot as plt
from src.fft_from_image.FFT import FFT
from src.eig_problem.EigenValueProblem import EigenValueProblem2D, EigenValueProblem1D
from src.io.DataReader import ParsingData
from src.modes.MagnetizationProfile import Profile2D

from src.eig_problem.InputParameter import InputParameter
from src.drawing.Plot import Plot

input_parameters = ParsingData('./tst/Parameter_for_TheImpact.yaml')
direction = 'x'
mode_number = 2
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

    eig_vec = EigenValueProblem2D(direction, input_parameters, 'Fe', 'Ni').calculate_eigen_vectors(bloch_vector=np.array([1,1]))
    np.savetxt('vec.vec', eig_vec.view(float))
    if show_plot:
        ax1 = plt.axes()
        Profile2D(50, eig_vec).generate_plot(mode_number)
    else:
        return 0


if __name__ == "__main__":
    do_program()
