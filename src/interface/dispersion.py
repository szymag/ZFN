import numpy as np
from src.fft_from_image.FFT import FFT
from src.eig_problem.EigenValueProblem import EigenValueProblem2D, EigenValueProblem1D
from src.io.DataReader import ParsingData
from src.drawing.Plot import Plot

input_parameters = ParsingData('./tst/Parameter_for_TheImpact.yaml')
direction = 'x'
number_of_dispersion_branch = 5
x_lim = None
y_lim = None
show_plot = True


def do_program():  # TODO: file type should represent containing data
    # TODO: Input file types hould be moved to DataReader
    if input_parameters.input_fft_file()[-3:] == 'png':
        fft = FFT().zwroc_tablice(input_parameters.input_fft_file())
        np.savetxt('tmp.fft', fft)
        input_parameters.set_new_value('tmp.vec', 'numerical_parameters', 'fft_file')
    elif input_parameters.input_fft_file()[-3:] == 'fft':
        fft = np.loadtxt(input_parameters.input_fft_file())
    elif input_parameters.input_fft_file()[-3:] == 'dys':
        return Plot(number_of_dispersion_branch, x_lim, y_lim).dispersion_relation(input_parameters.input_fft_file())
    else:
        raise ValueError

    eig_freq = EigenValueProblem2D(direction,
                                   input_parameters, 'Py', 'Co').calculate_dispersion_along_direction()
    np.savetxt(input_parameters.output_file('dispersion'), np.array(eig_freq))
    if show_plot:
        return Plot(number_of_dispersion_branch, x_lim, y_lim).dispersion_relation(eig_freq)
    else:
        return 0


def do_program_1D():
    eig_freq = EigenValueProblem1D(input_parameters, 'Co', 'Py').calculate_dispersion_along_direction()
    np.savetxt('dispersion_sokolovskyy_0.2T.tst', eig_freq)
    if show_plot:
        return Plot(number_of_dispersion_branch, x_lim, y_lim).dispersion_relation(input_parameters.input_fft_file())
    else:
        return 0


def do_program_oblique():
    eig_freq = EigenValueProblem1D('alongwires.yaml', 'Co', 'Py').oblique_dispersion()
    np.savetxt('11dispersion_sokolovskyy_0.2T.tst', eig_freq)
    Plot(5, x_lim, y_lim).dispersion_relation(input_parameters.output_file('dispersion'))


def do_program_map():
    eig_freq = EigenValueProblem2D(direction,
                                   input_parameters, 'Py', 'Co').calculate_dispersion_map()
    np.savetxt(input_parameters.output_file('dispersion'), eig_freq)

    if show_plot:
        return Plot(number_of_dispersion_branch, x_lim, y_lim).contour_plot(eig_freq, 20e9)
    else:
        return 0


if __name__ == "__main__":
    # do_program_oblique()
    do_program_map()
    #do_program_1D()
