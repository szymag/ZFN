import numpy as np
import matplotlib.pyplot as plt
from src.fft_from_image.FFT import FFT
from src.fft_from_image.Sequences import Phason, Fibonacci, Periodic
from src.io.DataReader import ParsingData
from src.modes.MagnetizationProfile import Profile1D
from src.interface.eigen_vector import do_program_1D
from src.interface.dispersion import do_program_idos
from src.drawing.Plot import Plot


input_parameters = ParsingData('./src/interface/Rychly.yaml')

repetition_seq = 5
fib_number = 11
green_stripes = 'Py'
gray_stripes = 'Co'
bloch_vec = [1e4, 0]

samples_count = 100
phasons = 10
sample_number = 0
sequence_type = 'F'


def define_name_phason(phasons_percentage, sample_number, structure_type, path=None):
    if path is not None:
        return path + structure_type + '_' + str(phasons_percentage) + '_' + str(sample_number)
    else:
        return structure_type + '_' + str(phasons_percentage) + '_' + str(sample_number)


def save_structure(phasons_percentage, samples_count, structure_type):
    for i in range(samples_count):
        f = Phason(structure_type, repetition_seq, fib_number, phasons_percentage)
        seq, phason_positions = f.sequence_shuffling(f.seq)

        seq = f.sequence(seq)
        fft_seq = FFT().fft1d(seq)
        np.savetxt(define_name_phason(phasons_percentage, i, structure_type) + '.fft', fft_seq)
        np.savetxt(define_name_phason(phasons_percentage, i, structure_type) + '.pos', phason_positions)


def save_eig_vector(file_name):
    fft_file = file_name
    input_parameters.set_new_value(fft_file + '.fft', 'numerical_parameters', 'fft_file')
    input_parameters.set_new_value(file_name + '.vec', 'numerical_parameters', 'output_file')
    print(input_parameters.output_file(''))
    do_program_1D(input_parameters, green_stripes, gray_stripes, bloch_vec)


def plot_modes(mode_number, file_name):
    plt.rc('xtick', labelsize='xx-large')
    plt.rc('ytick', labelsize='xx-large')
    ax1 = plt.axes()
    draw_structure(file_name, ax1)
    mode(mode_number, file_name).generate_plot(ax1, 0)


def draw_structure(file_name, axis):
    phasons = np.array(np.loadtxt(file_name + '.pos'), dtype=int)
    fib = Fibonacci(repetition_seq, fib_number)
    seq = fib.sequence_generator()
    Plot(1).draw_structure(axis, seq, phasons, 91)


def mode(mode_number, file_name):
    return Profile1D(mode_number, file_name + '.vec' , None, input_parameters)


def calculate_idos(file_name):
    input_parameters.set_new_value(file_name + '.fft', 'numerical_parameters', 'fft_file')
    input_parameters.set_new_value(file_name + '.dys', 'numerical_parameters', 'output_file')
    do_program_idos(input_parameters, green_stripes, gray_stripes, bloch_vec)


def plot_idos(phasons_percentage, start_point, end_point, structure_type, orginal_struct, path=None):
    plt.style.use('seaborn')
    plt.rc('text', usetex=False)
    plt.rc('xtick', labelsize='x-large')
    plt.rc('ytick', labelsize='x-large')
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel('Frequency (GHz)', fontsize='x-large')
    ax.set_ylabel('IDOS', fontsize='x-large')

    fib = np.loadtxt(orginal_struct)

    axins1 = ax.inset_axes([0.9, 0.9, 0.07, 0.07])
    idos_inset(ax, axins1, [19.5, 23], [270, 390], 0.05)

    axins3 = ax.inset_axes([0.9, 0.9, 0.05, 0.08])
    idos_inset(ax, axins3, [14.3, 16], [80, 100], 0.05)

    for i in range(start_point, end_point):
        a = np.loadtxt(define_name_phason(phasons_percentage, i, structure_type, path) + '.dys')
        for j in [axins1, axins3, ax]:
            Plot(1).idos(j, a, 'C0', alpha=0.02)

    for j in [ax, axins1, axins3]:
        Plot(1).idos(j, fib, 'C2', alpha=1)

    #ax.set_aspect(aspect=0.05)
    plt.tight_layout()


def idos_inset(axis, axis_inset, xlim, ylim, aspect):
    axis.indicate_inset_zoom(axis_inset)
    axis_inset.set_xlim(xlim)
    axis_inset.set_ylim(ylim)
    axis_inset.set_xticklabels('')
    axis_inset.set_yticklabels('')
    axis_inset.set_aspect(aspect=aspect)


def calculate_localization(mode_number, file_name, grid):
    mod = mode(mode_number, file_name).spatial_distribution_dynamic_magnetization(grid, mode_number)[1]
    mod =  mod / np.sum(abs(mod)) * grid
    return 1 / grid * np.sum(np.log(abs(mod)))


def calculate_fmr(mode_number, file_name, grid):
    mod_class = mode(mode_number, file_name)
    mod = mod_class.spatial_distribution_dynamic_magnetization(grid, mode_number)[1]
    return mod_class.fmr_intensity(mod)


def calculate_gap_statisitc(phasons_percentage, start_point, end_point, structure_type, range_to_look):
    gaps = np.zeros(end_point - start_point)
    for index, i in enumerate(range(start_point, end_point)):
        a = np.loadtxt(define_name_phason(phasons_percentage, i, 'F', path) + '.dys')
        gaps[index] = calculate_gap(a, range_to_look)
    return np.average(gaps), np.std(gaps)


def calculate_gap(data, range_to_look):
    return np.max(np.diff(data[range_to_look[0]:range_to_look[1]]))

if __name__ == "__main__":
    """
    calculate fft from disturbed structure
    """
    # save_structure(phasons, samples_count, 'F')

    """
    calculate modes
    """
    # file = define_name_phason(phasons, sample_number, sequence_type)
    # save_eig_vector('./f_coef_5*11')
    """
    Plot modes
    """
    # file = define_name_phason(phasons, 0, 'F')
    # for i in range(0, 9):
    #     plot_modes(i, file)
    #     plt.show()

    """
    calculate idos structure
    """
    # input_parameters.set_new_value('./f_coef_5*11.txt', 'numerical_parameters', 'fft_file')
    # input_parameters.set_new_value('F_5_11' + '.dys', 'numerical_parameters', 'output_file')
    # do_program_idos(input_parameters, green_stripes, gray_stripes, bloch_vec)
    #
    # for i in range(samples_count):
    #     file = define_name_phason(phasons, i, 'F')
    #     calculate_idos(file)
    """
    plotting idos structure
    """
    # plot_idos(phasons, 0, 99, 'F', './F_5_11.dys')
    # plt.show()
    # plt.savefig('idos.svg')
    """
    calculate localization factor
    """
    # a = np.zeros(100)
    # for i in range(100):
    #     # name = define_name_phason(phasons, 0, 'F')
    #     a[i] = calculate_localization(i, './f_coef_5*11', 1000)
    # plt.plot(a, np.arange(100))
    # print(a)
    # plt.show()
    """
    calculate fmr
    """
    a = np.zeros(100)
    for i in range(100):
        name = define_name_phason(phasons, 0, 'F')
        a[i] = calculate_fmr(i, name, 1000)
    plt.plot(np.arange(100), a)
    print(a)
    plt.show()
