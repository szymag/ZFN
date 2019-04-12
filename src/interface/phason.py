import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from src.fft_from_image.FFT import FFT
from src.eig_problem.EigenValueProblem import EigenValueProblem1D
from src.fft_from_image.Sequences import Phason, Fibonacci
from src.io.DataReader import ParsingData
from src.modes.MagnetizationProfile import Profile1D


input_parameters = ParsingData('./src/interface/Rychly.yaml')

repetition_seq = 10
fib_number = 14


def save_structure(phasons_percentage, samples_count):
    for i in range(samples_count):
        f = Phason('P', repetition_seq, fib_number, phasons_percentage)
        seq, phason_positions = f.sequence_shuffling(f.seq)

        seq = f.sequence(seq)
        fft_seq = FFT().fft1d(seq)
        np.savetxt('periodic_' + str(phasons_percentage) + '_' + str(i) + '.fft', fft_seq)
        np.savetxt('periodic_' + str(phasons_percentage) + '_' + str(i) + '.pos', phason_positions)


def save_eig_vector(phasons_percentage, sample_number):
    np.savetxt('phason_' + str(phasons_percentage) + '_' + str(sample_number) + '.vec',
               calculate_eig_vector(phasons_percentage, sample_number).view(float))


def calculate_eig_vector(phasons_percentage, sample_number):
    fft_file = './phason_' + str(phasons_percentage) + '_' + str(sample_number) + '.fft'
    print(fft_file)
    input_parameters.set_new_value(fft_file, 'numerical_parameters', 'fft_file')
    return EigenValueProblem1D(input_parameters, 'Py', 'Co').calculate_eigen_vectors(bloch_vector=[1e4, 0])


def mode(mode_number, phasons_percentage, sample_number, ax):
    vec_name = './phason_' + str(phasons_percentage) + '_' + str(sample_number)+'.vec'
    print(vec_name)
    mod = Profile1D(mode_number, vec_name,
                    None, input_parameters)
    return mod.generate_plot(ax, 0)


def draw_structure(phasons_percentage, sample_number, ax):
    phasons = np.array(np.loadtxt('./phason_' + str(phasons_percentage) + '_' + str(sample_number) + '.pos'), dtype=int) + 1
    fib = Fibonacci(repetition_seq, fib_number)
    seq = fib.sequence_generator()
    for index, el in enumerate(seq):
        if index in phasons:
            ax.add_patch(Rectangle((index * 91, 0), 91, 1, color='red', alpha=0.8, linewidth=0, edgecolor=None))
        elif index + 1 in phasons:
            ax.add_patch(Rectangle((index * 91, 0), 91, 1, color='red', alpha=0.4, linewidth=0, edgecolor=None))

        else:
            if el == 0:
                ax.add_patch(Rectangle((index*91, 0), 91, 1, color='green', alpha=0.5, linewidth=0, edgecolor=None))
            else:
                ax.add_patch(Rectangle((index*91, 0), 91, 1, color='gray', alpha=0.5, linewidth=0, edgecolor=None))


def plot_idos(ran):
    plt.style.use('seaborn')
    plt.rc('text', usetex=False)
    #plt.rc('font', family='serif')
    plt.rc('xtick',labelsize='large')
    plt.rc('ytick', labelsize='x-large')

    #plt.rcParams['axes.facecolor'] = '#F6FBFC'
    fig = plt.figure()
    ax = fig.add_subplot(111)
    fib = np.loadtxt('./idos_fib_10_14.dys')

    #average = np.zeros(len(fib))
    for i in range(ran[0], ran[1]):
        a = np.loadtxt('./idos_phason_0.1_' + str(i) + '.dys')
        ax.step(a / 1e9, np.arange(len(a)), alpha=0.02, color='C0')
    #ax.step(average / (ran[1] - ran[0]), np.arange(len(average)), color='C1')

    ax.step(fib / 1e9, np.arange(len(fib)), color='C2')
    ax.set_xlabel('Frequency (GHz)', fontsize='x-large')
    ax.set_ylabel('IDOS', fontsize='x-large')

    axins1 = ax.inset_axes([0.1, 0.5, 0.47, 0.47])
    for i in range(ran[0], ran[1]):
        a = np.loadtxt('./idos_phason_0.1_' + str(i) + '.dys')
        axins1.step(a / 1e9, np.arange(len(a)), alpha=0.02, color='C0')
    axins1.step(fib / 1e9, np.arange(len(fib)), color='C2')
    axins1.set_xlim(19.5, 23)
    axins1.set_ylim(270, 390)
    axins1.set_xticklabels('')
    axins1.set_yticklabels('')
    ax.indicate_inset_zoom(axins1)
    axins1.set_aspect(aspect=0.05)

    axins2 = ax.inset_axes([0.55, 0.05, 0.42, 0.3])
    for i in range(ran[0], ran[1]):
        a = np.loadtxt('./idos_phason_0.1_' + str(i) + '.dys')
        axins2.step(a / 1e9, np.arange(len(a)), alpha=0.02, color='C0')
    axins2.step(fib / 1e9, np.arange(len(fib)), color='C2')
    axins2.set_xlim(15, 20)
    axins2.set_ylim(130, 170)
    axins2.set_xticklabels('')
    axins2.set_yticklabels('')
    ax.indicate_inset_zoom(axins2)
    axins2.set_aspect(aspect=0.05)

    axins3 = ax.inset_axes([0.03, 0.15, 0.25, 0.28])
    for i in range(ran[0], ran[1]):
        a = np.loadtxt('./idos_phason_0.1_' + str(i) + '.dys')
        axins3.step(a / 1e9, np.arange(len(a)), alpha=0.02, color='C0')
    axins3.step(fib / 1e9, np.arange(len(fib)), color='C2')
    axins3.set_xlim(14.3, 16)
    axins3.set_ylim(80, 100)
    axins3.set_xticklabels('')
    axins3.set_yticklabels('')
    ax.indicate_inset_zoom(axins3)
    axins3.set_aspect(aspect=0.05)

    ax.set_aspect(aspect=0.05)
    plt.tight_layout()



if __name__ == "__main__":
    """
    calculate fft from disturbed structure
    """
    #save_structure(0.1, 100)

    """
    calculate modes
    """
    #save_eig_vector(0.1, 91)
    for i in range(50, 352):
        ax1 = plt.axes()
        draw_structure(0.1, 91, ax1)
        mode(i, 0.1, 91, ax1)
        plt.show()
        plt.clf()
        plt.close()
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
    # plot_idos((0, 99))
    # plt.savefig('idos.svg')
