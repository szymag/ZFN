import numpy as np
import matplotlib.pyplot as plt
from src.fft_from_image.FFT import FFT
from src.fft_from_image.Sequences import Phason, Fibonacci, Periodic, Random
from src.io.DataReader import ParsingData
from src.modes.MagnetizationProfile import Profile1D
from src.interface.eigen_vector import do_program_1D
from src.interface.dispersion import do_program_idos
from src.drawing.Plot import Plot
from src.utils.cProfiler import do_cprofile
from matplotlib.colors import LogNorm
from scipy.ndimage.filters import gaussian_filter
import matplotlib
from matplotlib.patches import Rectangle
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
    print('fft file: ' +str(fft_file), 'output: ' + str(input_parameters.output_file('')))
    do_program_1D(input_parameters, green_stripes, gray_stripes, bloch_vec)


def plot_modes(mode_number, file_name, structure_type, freq):
    plt.rc('xtick', labelsize='small')
    plt.rc('ytick', labelsize='small')
    plt.style.use('seaborn')
    plt.rc('text', usetex=False)
    fig, axs = plt.subplots(4,figsize=(5, 5))
    print(axs)
    size = 91*377/4/1000
    for index, i in enumerate(axs):
        print(size*index, size*(index+1))
        i.set_xlim(size*index, size*(index+1))
        draw_structure(file_name, i, structure_type)
        mode(mode_number, file_name).generate_plot(i, 0)
    axs[0].set_title('mode: ' + str(mode_number) + '\n'
                     +'frequency: ' + str(round(freq/1e9,3)) + ' GHz')


def draw_structure(file_name, axis, structure_type):
    print(file_name)
    phasons = np.array(np.loadtxt(file_name + '.pos'), dtype=int)
    if structure_type == 'F':
        struct = Fibonacci(repetition_seq, fib_number)
    elif structure_type == "P":
        struct = Periodic(repetition_seq, fib_number)
    elif structure_type == "R":
        struct = Random(repetition_seq, Fibonacci(repetition_seq, fib_number).fib_number(),
                         Fibonacci(repetition_seq, fib_number- 2).fib_number())
    seq = struct.sequence_generator()
    Plot(1).draw_structure(axis, seq, phasons, 91/1000)


def mode(mode_number, file_name):
    return Profile1D(mode_number, file_name + '.vec', None, input_parameters)


def calculate_idos(file_name):
    input_parameters.set_new_value(file_name + '.fft', 'numerical_parameters', 'fft_file')
    input_parameters.set_new_value(file_name + '.dys', 'numerical_parameters', 'output_file')
    do_program_idos(input_parameters, green_stripes, gray_stripes, bloch_vec)


def plot_idos(phasons_percentage, start_point, end_point, structure_type, orginal_struct, path=None):
    plt.style.use('seaborn-colorblind')
    plt.rc('text', usetex=False)
    # plt.rc('xtick', labelsize='x-large')
    # plt.rc('ytick', labelsize='x-large')
    # fig = plt.figure()
    fig = plt.figure(figsize=(3.3,5))
    ax = fig.add_subplot(111)
    ax.set_xlabel('Frequency (GHz)')
    ax.set_ylabel('IDOS')

    fib = np.loadtxt(orginal_struct)

    # axins1 = ax.inset_axes([0.2, 0.7, 0.4, 0.3])
    # idos_inset(axins1, [19.5, 23], [270, 300], 0.2)
    #
    # axins3 = ax.inset_axes([0.5, 0.01, 0.4, 0.4])
    # idos_inset(axins3, [14.5, 16], [75, 115], 0.02)

    # for i in range(start_point, end_point):
    #     a = np.loadtxt(define_name_phason(phasons_percentage, i, structure_type, path) + '.dys')
    #     for j in [axins1, axins3, ax]:
    #         Plot(1).idos(j, a, '#0A539E', alpha=1)
    a = np.loadtxt(define_name_phason(phasons_percentage, i, structure_type, path) + '.dys')
    Plot(1).idos(ax, a, '#0A539E', alpha=1)
    # for j in [axins1, axins3, ax]:
    #     Plot(1).idos(j, fib, 'C2', alpha=1)

    # ax.set_aspect(aspect=0.05)
    plt.tight_layout()

def color_generator():
    yield 'C4'
    yield 'C1'
    yield 'C2'
    yield 'C3'
    while True:
        yield 'gray'

def add_gaps_idos(axis):
    x_axis = np.loadtxt('Fib_14.dys')/1e9
    x_axis = x_axis[30:]
    gaps = x_axis[np.argsort(np.diff(x_axis))[-11:]][::-1]
    gaps_width = np.sort(np.diff(x_axis))[-11:][::-1]
    color = color_generator()
    print(gaps_width)
    for g, g_w in zip(gaps, gaps_width):
        axis.add_patch(Rectangle((g, -4), g_w, 440.1,alpha=0.7, color=next(color)))


def idos_inset(axis_inset, xlim, ylim, aspect):
    axis_inset.set_xlim(xlim)
    axis_inset.set_ylim(ylim)
    axis_inset.set_xticklabels('')
    axis_inset.set_yticklabels('')
    axis_inset.set_aspect(aspect=aspect)


def calculate_localization(mode, grid):
    mod = mode / np.sum(abs(mode)) * grid
    return -1 / grid * np.sum(abs(mod)*np.log(abs(mod)))


def plot_localization(file_name, grid, title, type=None):
    plt.style.use('seaborn')
    plt.rc('text', usetex=False)
    plt.rc('xtick', labelsize='xx-large')
    plt.rc('ytick', labelsize='xx-large')
    a = np.zeros(400)
    mod = mode(0, file_name)
    for index in range(len(a)):
        a[index] = calculate_localization(mod.spatial_distribution_dynamic_magnetization(grid, index)[1], grid) # calculation
        np.savetxt(title + '.txt', a)
    # a = np.loadtxt(title + ".txt")
    currentAxis = plt.gca()
    add_gaps(currentAxis, file_name)
    x_label = np.loadtxt(file_name + '.dys')[0:len(a)]
    if type == "ref":
        plt.scatter(x_label/1e9, a, s=15, alpha=1, c="C1", linewidth=None, edgecolors=None)
    else:
        plt.scatter(x_label/1e9, a, s=15, alpha=0.2, c="C0", linewidth=None, edgecolors=None)

    plt.ylabel('$\lambda$', fontsize='xx-large')
    plt.xlabel('Frequency (GHz)', fontsize='xx-large')
    plt.ylim([-4, 0.1])
    plt.xlim([10,25])
    # plt.title(title)



def add_gaps(axis, file_name):
    x_axis = np.loadtxt('Fib_14.dys')/1e9
    x_axis = x_axis[20:]
    gaps = x_axis[np.argsort(np.diff(x_axis))[-11:]]
    gaps_width = np.sort(np.diff(x_axis))[-11:]
    for g, g_w in zip(gaps, gaps_width):
        axis.add_patch(Rectangle((g, -4), g_w, 4.1,alpha=0.1 , color="C2"))


def plot_fmr(file_name, grid, title):
    plt.style.use('seaborn')
    plt.rc('text', usetex=False)
    plt.rc('xtick', labelsize='xx-large')
    plt.rc('ytick', labelsize='xx-large')
    a = np.zeros(400)

    mod = mode(0, file_name)
    for index in range(len(a)):
        a[index] = calculate_fmr(mod.spatial_distribution_dynamic_magnetization(grid, index)[1])

    x_label = np.loadtxt(file_name + '.dys')[0:len(a)]
    a[0] = 1000
    plt.scatter(x_label, a, s=25)
    plt.ylabel('Intensity', fontsize='xx-large')
    plt.xlabel('Frequency (GHz)', fontsize='xx-large')
    print(np.argsort(a)[-5:][::-1])
    plt.title(title)
    plt.savefig(title + '.svg')


def calculate_fmr(mode):
    return np.abs(np.trapz(mode)) ** 2 / np.trapz(np.abs(mode) ** 2)


def plot_gap_width(phasons, frequency_ranges, title):
    plt.style.use('seaborn')
    plt.rc('text', usetex=False)
    plt.rc('xtick', labelsize='xx-large')
    plt.rc('ytick', labelsize='xx-large')
    results = np.zeros((4, len(phasons)))
    error = np.zeros((4, len(phasons)))
    color = color_generator()
    for num, freq in enumerate(frequency_ranges):
        col = next(color)
        for index, i in enumerate(phasons):
            results[num][index], error[num][index] = calculate_gap_statisitc(i, 0, 100, 'F', freq)
        plt.errorbar(phasons, results[num], error[num], fmt='--o', label='gap at: ' + str((freq[0] + freq[1])/2/1e9) + ' GHz', color=col)
    plt.legend()
    plt.xlabel("Phasons", fontsize='xx-large')
    plt.ylabel("Gap width (GHz)", fontsize='xx-large')
    plt.savefig(title + '.svg')


def calculate_gap_statisitc(phasons_percentage, start_point, end_point, structure_type, range_to_look):
    gaps = np.zeros(end_point - start_point)
    for index, i in enumerate(range(start_point, end_point)):
        a = np.loadtxt(define_name_phason(phasons_percentage, i, structure_type, './JEMS2019/fib/') + '.dys')
        gaps[index] = calculate_gap(a, range_to_look)
    return np.average(gaps)/1e9, np.std(gaps)/1e9


def calculate_gap(data, range_to_look):
    range_to_look = [np.argmin(abs(data - i)) for i in range_to_look]
    return np.max(np.diff(data[range_to_look[0]:range_to_look[1]]))

def idos_histogram(phasons_percentage, structure_type, path):
    plt.style.use('seaborn-colorblind')
    plt.rc('text', usetex=False)
    matplotlib.rcParams['font.family'] = "Liberation Sans"

    # file_x = [np.loadtxt(define_name_phason(phasons_percentage, 2, structure_type, path) + '.dys')/1e9  for i in range(100)]
    file_x = [np.loadtxt('./Fib_14.dys')/1e9  for i in range(100)]
    file_y = [np.arange(0, 400)  for i in range(100)]

    x = np.concatenate(file_x, axis=None)
    y = np.concatenate(file_y, axis=None)

    h, x, y= np.histogram2d(x, y, bins=(len(x)/100,len(y)/100))
    X, Y = np.meshgrid(x, y)
    h = h.T/np.max(h)
    from src.interface.bragg import dispersion
    k = np.linspace(1e1,0.4e8, 400)
    lattice_const = 34307e-9
    fig = plt.figure(figsize=(3.3,4.4))
    ax2 = fig.add_subplot(211)
    ax = ax2.twiny()
    # ax2.set_ylim([0, 0.4e8/2/np.pi*lattice_const])
    ax2.plot(k/2/np.pi*lattice_const, dispersion(k), color='C2', alpha=0)
    ax.plot(np.arange(0,400), dispersion(k), color='C2')
    ax2.tick_params(axis='x', labelcolor='C2')
    ax.tick_params(axis='x', labelcolor='C0')
    ax2.set_xlabel('Wavenumber (2$\pi$/a)', color="C2")
    # Plot(1).idos(ax, np.loadtxt('./Fib_14.dys'), 'C2', alpha=0.5)

    a = np.arange(400)
    # x_label = np.loadtxt('./Fib_14.dys')


    # plt.scatter(x_label/1e9, a, s=0.7, alpha=1, c="C2", linewidth=None, edgecolors=None)
    ax2.set_ylabel('Frequency (GHz)')
    ax.set_xlabel('IDOS', color='C0')
    # plt.ylim([-4, 0.1])
    # plt.xlim([10,22.5])

    ax.pcolormesh(Y, X, gaussian_filter(h, sigma=1, order=0), cmap="Blues", norm=LogNorm(vmin=0.001, vmax=1))
    # ax.grid(True, linewidth=0.5)


    ax.set_ylim([9.98, 24])


    # ax.plot(np.arange(400), dispersion(k), color='C2')
    ax.set_ylim([9.98, 24])
    ax.set_xlim([-15, 410])
    plt.tight_layout()
    # add_gaps_vertical(ax)

    ax1 = fig.add_subplot(212)
    ax1.set_ylabel('Intenisty (arb. units)')
    ax1.set_xlabel('Wavenumber (2$\pi$/a)', color="C2")
    ax1.tick_params(axis='x', labelcolor='C2')
    full_peak = abs(np.loadtxt('./f_coef_10*14.fft').view(complex))
    full_peak = full_peak[len(full_peak)//2:]/max(full_peak)

    coef = np.transpose(np.loadtxt('/media/szymag/Dane/ZFN/JEMS2019/fib/F_'+ str(phasons_percentage) +'_0.fft'))
    coef = coef[0] + coef[1]*1j
    struct = abs(np.fft.ifft(coef))
    coef = coef[len(coef)//2:]/max(coef)
    coef = abs(coef[:])

    x_coef = 2*np.pi*np.arange(0, len(coef))/lattice_const
    ax1.set_xlim(-0.35e7/2/np.pi*lattice_const, 0.8e8/2/np.pi*lattice_const)
    ax1.plot(np.arange(0, len(coef)), full_peak, color='C2')
    # ax1.plot(np.arange(0, len(coef)), coef[:len(coef)], color='C2')
    plt.tight_layout()

def add_gaps(axis):
    x_axis = np.loadtxt('Fib_14.dys')/1e9
    x_axis = x_axis[20:]
    gaps = x_axis[np.argsort(np.diff(x_axis))[-11:]][::-1]
    gaps_width = np.sort(np.diff(x_axis))[-11:][::-1]
    for g, g_w in zip(gaps, gaps_width):
        axis.add_patch(Rectangle((g, -4), g_w, 410,alpha=0.3 , color='gray', ec=None))

def add_gaps_vertical(axis):
    y_axis = np.loadtxt('Fib_14.dys')/1e9
    y_axis = y_axis[20:]
    gaps = y_axis[np.argsort(np.diff(y_axis))[-11:]][::-1]
    gaps_width = np.sort(np.diff(y_axis))[-11:][::-1]
    for g, g_w in zip(gaps, gaps_width):
        axis.add_patch(Rectangle((-4, g), 410,  g_w, alpha=0.3 , color='gray', ec=None))


if __name__ == "__main__":
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
    """
    calculate fft from disturbed structure
    """
    # for i in [1, 2]:
    #    save_structure(i, samples_count, sequence_type)
    # save_structure(144, samples_count, sequence_type)

    """
    calculate modes
    """
    # for i in [1, 2]:
    #     for j in range(5):
    #         file = define_name_phason(i, j, sequence_type)
    #         save_eig_vector(file)
    # save_eig_vector('./f_coef_10*14')

    """
    Plot modes
    """
    # file = define_name_phason(25, 3, 'F', '/media/szymag/Dane/ZFN/JEMS2019/fib/')
    # for i in [300, 302, 288, 291]:
    #      plot_modes(i, file)
    #      plt.show()

    # ax1 = plt.axes()
    # file_name = '/media/szymag/Dane/ZFN/JEMS2019/fib/F_15_3'
    # idos = np.loadtxt('/media/szymag/Dane/ZFN/JEMS2019/fib/F_15_3.dys')
    # # draw_structure(file_name, ax1)
    # for i in range(30, 50):
    #     plot_modes(i, file_name, 'F', idos[i])
    # # ax1.set_xlim(10000)
    # # plt.savefig('emty.svg')
    #     plt.savefig('./15p_3/mode_' + str(i) +'.svg')
    #     plt.close()
    #     plt.clf()
    """
    calculate idos structure
    """

    # for phas in [15]:
    #     idos_histogram(phas, 'F', '/media/szymag/Dane/ZFN/JEMS2019/fib/')
    #     plt.xlabel('Frequency (GHz)')
    #     plt.ylabel('IDOS')
    #     plt.savefig('idos_' + str(phas) + '.svg')
    #     plt.clf()
    #     plt.cla()
    # input_parameters.set_new_value('./f_coef_10*14.txt', 'numerical_parameters', 'fft_file')
    # input_parameters.set_new_value('Fib_14' + '.dys', 'numerical_parameters', 'output_file')
    # do_program_idos(input_parameters, green_stripes, gray_stripes, bloch_vec)
    # #
    # for phas in [5, 35, 65]:
    #
    #     for i in range(samples_count):
    #         file = define_name_phason(phas, i, 'F', './JEMS2019/fib/')
    #         print(file)
    #         calculate_idos(file)
    # for i in range(samples_count):
    #     file = define_name_phason(144, i, 'R')
    #     calculate_idos(file)
    """
    plotting idos structure
    """
    # for i in [5, 15, 25, 35, 45, 55, 65, 80, 100, 120, 140]:
    for i in [45]:
        idos_histogram(i, 'F', '/media/szymag/Dane/ZFN/JEMS2019/fib/')
    # plt.show()

    plt.savefig('bragg.png', dpi=400)

    # mode_number = 50
    # idos = np.loadtxt('/media/szymag/Dane/ZFN/JEMS2019/fib/F_15_3.dys')
    # for mode_number in range(400):
    #
    #     plt.style.use('seaborn-colorblind')
    #     plt.rc('text', usetex=False)
    #     fig = plt.figure(figsize=(10,7))
    #     gs = fig.add_gridspec(4,2)
    #     ax1 = fig.add_subplot(gs[:, 0])
    #     ax2 = fig.add_subplot(gs[0, 1])
    #     ax3 = fig.add_subplot(gs[1, 1])
    #     ax4 = fig.add_subplot(gs[2, 1])
    #     ax5 = fig.add_subplot(gs[3, 1])
    #     # ax = fig.add_subplot(121)
    #     ax1.set_xlabel('Frequency (GHz)')
    #     ax1.set_ylabel('IDOS')
    #
    #     a = np.loadtxt(define_name_phason(15, 3, 'F', '/media/szymag/Dane/ZFN/JEMS2019/fib/') + '.dys')
    #     Plot(1).idos(ax1, a, '#0A539E', alpha=1)
    #
    #
    #
    #     axs = [ax2, ax3, ax4, ax5]
    #     size = 91*377/4/1000
    #
    #     ax1.scatter(idos[mode_number]/1e9, mode_number, s=100, color='r')
    #     for index, i in enumerate(axs):
    #         i.set_xlim(size*index, size*(index+1))
    #         draw_structure('/media/szymag/Dane/ZFN/JEMS2019/fib/F_15_3', i, 'F')
    #         mode(mode_number, '/media/szymag/Dane/ZFN/JEMS2019/fib/F_15_3').generate_plot(i, 0)
    #     axs[0].set_title('mode: ' + str(mode_number) + '\n'
    #                      +'frequency: ' + str(round(idos[mode_number]/1e9,3)) + ' GHz')
    #     plt.tight_layout()
    #     plt.savefig(str(mode_number) + '.png')
    #     plt.close(fig)
        # plt.clf()
    # plt.show()

    #     plt.savefig('idos_' + str(i) + '.svg')
    #     plt.clf()
    #     plt.close()
    # plot_idos(144, 0, 100, 'R', './Fib_14.dys')
    # plt.savefig('idos_' + 'random' + '.svg')
    """
    plot localization factor
    """
    # import argparse
    # parser = argparse.ArgumentParser()
    # parser.add_argument("count", type=int)
    # phas = parser.parse_args().count
    # print(phas)
    # for i in range(5):
    #    name = define_name_phason(phas, i, 'F')
    #    print(name)
    #    plot_localization(name, 7000, 'lambda, ' + str(phas) + ' phasons ' + str(i) + ' series')
    #    # plt.clf()
    # plot_localization('./Fib_14', 5000, 'lambda, 0 phason', 'ref')
    # plt.savefig('lambda_' + str(phas) + '.png', dpi=200)


    """
    calculate fmr
    """
    # for phas in [5, 15, 25, 35, 45, 55, 65, 80, 100, 120, 140]:
    #     for i in range(5):
    #         name = define_name_phason(phas, i, 'P')
    #         plot_fmr(name, 7000, 'fmr ' + str(phas) + ' phasons, '+ str(i) +  ' series')
    #     plt.clf()
        #plot_fmr('./Fib_14', 7000, '0 phasons')
    """
    gaps_width
    """
    a = [5, 15, 25, 35, 45, 55, 65, 80, 100, 120, 140]
    b = [(16.05e9, 19.4e9), (20.1e9, 21.7e9), (15e9, 15.5e9), (13.8e9, 14.2e9)]
    plot_gap_width(a, b, 'localization')
