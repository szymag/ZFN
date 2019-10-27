import numpy as np
import matplotlib.pyplot as plt
from src.fft_from_image.FFT import FFT
from src.fft_from_image.Sequences import Phason, Fibonacci, Periodic
from src.io.DataReader import ParsingData
from src.modes.MagnetizationProfile import Profile1D
from src.interface.eigen_vector import do_program_1D
from src.interface.dispersion import do_program_idos
from src.drawing.Plot import Plot
from src.utils.cProfiler import do_cprofile

input_parameters = ParsingData('./src/interface/Rychly.yaml')

repetition_seq = 10
fib_number = 14
green_stripes = 'Py'
gray_stripes = 'Co'
bloch_vec = [1e4, 0]

samples_count = 100
phasons = 5
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
    return Profile1D(mode_number, file_name + '.vec', None, input_parameters)


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

    axins1 = ax.inset_axes([0.2, 0.7, 0.4, 0.3])
    idos_inset(axins1, [19.5, 23], [270, 300], 0.2)

    axins3 = ax.inset_axes([0.5, 0.01, 0.4, 0.4])
    idos_inset(axins3, [14.5, 16], [75, 115], 0.02)

    for i in range(start_point, end_point):
        a = np.loadtxt(define_name_phason(phasons_percentage, i, structure_type, path) + '.dys')
        for j in [axins1, axins3, ax]:
            Plot(1).idos(j, a, 'C0', alpha=0.02)

    for j in [axins1, axins3, ax]:
        Plot(1).idos(j, fib, 'C2', alpha=1)

    ax.set_aspect(aspect=0.05)
    plt.tight_layout()


def idos_inset(axis_inset, xlim, ylim, aspect):
    axis_inset.set_xlim(xlim)
    axis_inset.set_ylim(ylim)
    axis_inset.set_xticklabels('')
    axis_inset.set_yticklabels('')
    axis_inset.set_aspect(aspect=aspect)


def calculate_localization(mode, grid):
    mod = mode / np.sum(abs(mode)) * grid
    return 1 / grid * np.sum(np.log(abs(mod)))


def plot_localization(file_name, grid, title):
    plt.style.use('seaborn')
    plt.rc('text', usetex=False)
    plt.rc('xtick', labelsize='x-large')
    plt.rc('ytick', labelsize='x-large')
    a = np.zeros(400)

    mod = mode(0, file_name)
    for index in range(len(a)):
        a[index] = calculate_localization(mod.spatial_distribution_dynamic_magnetization(grid, index)[1], grid)
        print(index)
    x_label = np.loadtxt(file_name + '.dys')[0:len(a)]
    plt.scatter(x_label, a, s=25)
    plt.ylabel('$\lambda$', fontsize='x-large')
    plt.xlabel('Frequency (GHz)', fontsize='x-large')
    plt.ylim([-8, 0.1])
    plt.title(title)
    plt.savefig(title + '.svg')


def plot_fmr(file_name, grid, title):
    plt.style.use('seaborn')
    plt.rc('text', usetex=False)
    plt.rc('xtick', labelsize='x-large')
    plt.rc('ytick', labelsize='x-large')
    a = np.zeros(400)

    mod = mode(0, file_name)
    for index in range(len(a)):
        a[index] = calculate_fmr(mod.spatial_distribution_dynamic_magnetization(grid, index)[1])
    x_label = np.loadtxt(file_name + '.dys')[0:len(a)]
    # a[0] = 1000
    plt.scatter(x_label, a, s=25)
    plt.ylabel('Intensity', fontsize='x-large')
    plt.xlabel('Frequency (GHz)', fontsize='x-large')
    print(np.argsort(a)[-5:][::-1])
    plt.title(title)
    plt.savefig(title + '.svg')


def calculate_fmr(mode):
    return np.abs(np.trapz(mode)) ** 2 / np.trapz(np.abs(mode) ** 2)


def plot_gap_width(phasons, frequency_ranges, title):
    plt.style.use('seaborn')
    plt.rc('text', usetex=False)
    plt.rc('xtick', labelsize='x-large')
    plt.rc('ytick', labelsize='x-large')
    results = np.zeros((4, len(phasons)))
    error = np.zeros((4, len(phasons)))
    for num, freq in enumerate(frequency_ranges):
        for index, i in enumerate(phasons):
            results[num][index], error[num][index] = calculate_gap_statisitc(i, 0, 100, 'F', freq)
        plt.errorbar(phasons, results[num], error[num], fmt='--o', label='gap at: ' + str((freq[0] + freq[1])/2/1e9) + ' GHz')
    plt.legend()
    plt.xlabel("Phasons", fontsize='x-large')
    plt.ylabel("Gap width (GHz)", fontsize='x-large')
    plt.savefig(title + '.svg')


def calculate_gap_statisitc(phasons_percentage, start_point, end_point, structure_type, range_to_look):
    gaps = np.zeros(end_point - start_point)
    for index, i in enumerate(range(start_point, end_point)):
        a = np.loadtxt(define_name_phason(phasons_percentage, i, structure_type) + '.dys')
        gaps[index] = calculate_gap(a, range_to_look)
    return np.average(gaps)/1e9, np.std(gaps)/1e9


def calculate_gap(data, range_to_look):
    range_to_look = [np.argmin(abs(data - i)) for i in range_to_look]
    return np.max(np.diff(data[range_to_look[0]:range_to_look[1]]))


if __name__ == "__main__":
    """
    calculate fft from disturbed structure
    """
    # for i in [5, 15, 25, 35, 45, 55, 65, 80, 100, 120, 140]:
    #    save_structure(i, samples_count, 'F')

    """
    calculate modes
    """
    # for i in [5, 15, 25, 35, 45, 55, 65, 80, 100, 120, 140]:
    #     for j in range(10):
    #         file = define_name_phason(phasons, sample_number, sequence_type)
    #         file = define_name_phason(i, j, 'F')
    #         save_eig_vector(file)
    # save_eig_vector('./f_coef_10*14')
    """
    Plot modes
    """
    # file = define_name_phason(25, 3, 'F')
    # for i in [300, 302, 288, 291]:
    #      plot_modes(i, file)
    #      plt.show()

    """
    calculate idos structure
    """
    # input_parameters.set_new_value('./f_coef_10*14.txt', 'numerical_parameters', 'fft_file')
    # input_parameters.set_new_value('Fib_14' + '.dys', 'numerical_parameters', 'output_file')
    # do_program_idos(input_parameters, green_stripes, gray_stripes, bloch_vec)
    # #
    # for phas in [5, 15, 25, 35, 45, 55, 65, 80, 100, 120, 140]:
    #     for i in range(samples_count):
    #         file = define_name_phason(phas, i, 'F')
    #         calculate_idos(file)
    """
    plotting idos structure
    """
    # for i in [5, 15, 25, 35, 45, 55, 65, 80, 100, 120, 140]:
    #     plot_idos(i, 0, 100, 'F', './Fib_14.dys')
    # # plt.show()
    #     plt.savefig('idos_' + str(i) + '.svg')
    #     plt.clf()
    #     plt.close()
    """
    plot localization factor
    """
    # for i in range(10):
    #    name = define_name_phason(140, i, 'F')
    #    plot_localization(name, 7000, 'lambda, 140 phason ' + str(i) + ' series')

    #plot_localization('./Fib_14', 5000, 'lambda, 0 phason')

    """
    calculate fmr
    """
    # for i in range(3, 4):
    #     name = define_name_phason(25, i, 'F')
    #     plot_fmr(name, 7000, '5 phasons, '+ str(i) +  ' series')
    # plot_fmr('./Fib_14', 7000, '0 phasons')
    """
    calculate localization
    """
    # a = [5, 15, 25, 35, 45, 55, 65, 80, 100, 120, 140]
    # b = [(13.8e9, 14.2e9), (15e9, 15.5e9), (16.05e9, 19.4e9), (20.1e9, 21.7e9)]
    # plot_gap_width(a, b, 'localization')
