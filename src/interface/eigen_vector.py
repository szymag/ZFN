import numpy as np
import matplotlib.pyplot as plt
from src.fft_from_image.FFT import FFT
from src.eig_problem.EigenValueProblem import EigenValueProblem2D, EigenValueProblem1D
from src.io.DataReader import ParsingData
from src.modes.MagnetizationProfile import Profile2D, Profile1D


from src.drawing.Plot import Plot

input_parameters = ParsingData('./src/interface/Rychly.yaml')
direction = 'x'
mode_number = 2
x_lim = None
y_lim = None
show_plot = True
# TODO: fft_data is wrong name. This argument is broader


def do_program():
    if input_parameters.fft_data()[-3:] == 'png':
        fft = FFT().zwroc_tablice(input_parameters.fft_data())
        np.savetxt('tmp.fft', fft)
        input_parameters.set_new_value('tmp.vec', 'numerical_parameters', 'fft_file')
    elif input_parameters.fft_data()[-3:] == 'fft':
        pass
    elif input_parameters.fft_data()[-3:] == 'vec':
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


def do_program_1D(input_param, mat_1, mat_2, bloch_vec):
    eig_vec = EigenValueProblem1D(input_param, mat_1, mat_2).calculate_eigen_vectors(bloch_vector=bloch_vec)
    np.savetxt(input_param.output_file('vectors'), eig_vec.view(float))


if __name__ == "__main__":
    pass
