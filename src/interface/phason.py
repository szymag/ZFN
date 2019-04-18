import numpy as np
import matplotlib.pyplot as plt
from src.fft_from_image.FFT import FFT
from src.fft_from_image.Sequences import Phason, Fibonacci
from src.io.DataReader import ParsingData
from src.modes.MagnetizationProfile import Profile1D
from src.interface.eigen_vector import do_program_1D
from src.drawing.Plot import Plot


input_parameters = ParsingData('./src/interface/Rychly.yaml')

repetition_seq = 10
fib_number = 14
green_stripes = 'Py'
gray_stripes = 'Co'
bloch_vec = [1e4, 0]


def save_structure(phasons_percentage, samples_count):
    for i in range(samples_count):
        f = Phason('P', repetition_seq, fib_number, phasons_percentage)
        seq, phason_positions = f.sequence_shuffling(f.seq)

        seq = f.sequence(seq)
        fft_seq = FFT().fft1d(seq)
        np.savetxt(define_name_phason(phasons_percentage, i) + '.fft', fft_seq)
        np.savetxt(define_name_phason(phasons_percentage, i) + '.pos', phason_positions)


def save_eig_vector(phasons_percentage, sample_number):
    fft_file = define_name_phason(phasons_percentage, sample_number) + '.fft'
    input_parameters.set_new_value(fft_file, 'numerical_parameters', 'fft_file')

    np.savetxt(define_name_phason(phasons_percentage, sample_number) + '.vec',
               do_program_1D(input_parameters, green_stripes, gray_stripes, bloch_vec))


def define_name_phason(phasons_percentage, sample_number, path=None):
    if path is not None:
        return path + 'phason_' + str(phasons_percentage) + '_' + str(sample_number)
    else:
        return 'phason_' + str(phasons_percentage) + '_' + str(sample_number)


def plot_modes():
    pass


def mode(mode_number, phasons_percentage, sample_number, ax):
    vec_name = define_name_phason(phasons_percentage, sample_number) + '.vec'
    print(vec_name)
    mod = Profile1D(mode_number, vec_name, None, input_parameters)
    return mod.generate_plot(ax, 0)


def draw_structure(phasons_percentage, sample_number, ax):
    phasons = np.array(np.loadtxt(define_name_phason(phasons_percentage, sample_number) + '.pos'), dtype=int) + 1
    fib = Fibonacci(repetition_seq, fib_number)
    seq = fib.sequence_generator()
    Plot(1).draw_structure(ax, seq, phasons, 91)


def plot_idos(phasons_percentage, start_point, end_point, path=None):
    plt.style.use('seaborn')
    plt.rc('text', usetex=False)
    plt.rc('xtick',labelsize='x-large')
    plt.rc('ytick', labelsize='x-large')
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel('Frequency (GHz)', fontsize='x-large')
    ax.set_ylabel('IDOS', fontsize='x-large')

    fib = np.loadtxt('./idos_fib_10_14.dys')

    axins1 = ax.inset_axes([0.1, 0.5, 0.47, 0.47])
    idos_inset(ax, axins1, [19.5, 23], [270, 390], 0.05)

    axins3 = ax.inset_axes([0.03, 0.15, 0.25, 0.28])
    idos_inset(ax, axins3, [14.3, 16], [80, 100], 0.05)

    for i in range(start_point, end_point):
        a = np.loadtxt(define_name_phason(phasons_percentage, i, path) + '.dys')
        for j in [axins1, axins3, ax]:
            Plot(1).idos(j, a, 'C0', alpha=0.02)

    for j in [ax, axins1, axins3]:
        Plot(1).idos(j, fib, 'C2', alpha=1)

    ax.set_aspect(aspect=0.05)
    plt.tight_layout()


def idos_inset(axis, axis_inset, xlim, ylim, aspect):
    axis.indicate_inset_zoom(axis_inset)
    axis_inset.set_xlim(xlim)
    axis_inset.set_ylim(ylim)
    axis_inset.set_xticklabels('')
    axis_inset.set_yticklabels('')
    axis_inset.set_aspect(aspect=aspect)


if __name__ == "__main__":
    """
    calculate fft from disturbed structure
    """
    #save_structure(0.1, 100)

    """
    calculate modes
    """
    #save_eig_vector(0.1, 91)
    # for i in range(50, 352):
    #     ax1 = plt.axes()
    #     draw_structure(0.1, 91, ax1)
    #     mode(i, 0.1, 91, ax1)
    #     plt.show()
    #     plt.clf()
    #     plt.close()
    # plt.rc('xtick', labelsize='xx-large')
    # plt.rc('ytick', labelsize='xx-large')
    # ax1 = plt.axes()
    # mode_num = 140
    # draw_structure(0.1, 91, ax1)
    # mode(mode_num, 0.1, 91, ax1)
    # ax1.set_xlim(15, 17.5)
    # ax1.set_aspect(aspect=0.4)
    # plt.savefig('mode ' + str(mode_num) + 'f.svg')
    # plt.clf()
    # plt.close()
    """
    plotting idos structure
    """
    plot_idos(0.1, 1, 99, './phason - JEMS2019/idos/')
    plt.show()
    # plt.savefig('idos.svg')
