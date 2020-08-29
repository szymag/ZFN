import numpy as np
import matplotlib.pyplot as plt
from src.fft_from_image.FFT import FFT
from src.fft_from_image.Sequences import Phason, Fibonacci
from src.io.DataReader import ParsingData
from src.modes.MagnetizationProfile import Profile1D
from src.interface.eigen_vector import do_program_1D
from src.interface.dispersion import do_program_idos
from src.eig_problem.EigenValueProblem import EigenValueProblem1D
from src.drawing.Plot import Plot
import matplotlib.patches as patches

import matplotlib
plt.style.use('seaborn-colorblind')
plt.rc('text', usetex=False)
matplotlib.rcParams['font.family'] = "Liberation Sans"

def generate_sequences(file_count=50, rec_vector_count=150):
    step = int(np.ceil(rec_vector_count / file_count))
    tab = []
    for i in range(file_count):
        arr = np.zeros(rec_vector_count)
        arr[:i*step] += 1
        tab.append(arr)
    return tab


def calculate_fft(sequence_list):
    for index, i in enumerate(sequence_list):
        # print(i)
        tab =  np.fft.fftshift(np.fft.fft(i))  / len(i)
        tab = np.stack((tab.real, tab.imag), axis=-1)
        np.savetxt('./src/interface/fft_topo/fft_' + str(index) + '.fft', tab)
    print('Files saved to ./fft_topo directory')


def draw_elementary_cell(fft_file):
    print(fft_file.shape)
    plt.plot(np.abs(fft_file))
    plt.plot(np.abs(np.fft.fft(fft_file)))
    plt.show()


def calculate_dispersion(ff1, ff2):
    for i in [ff1, ff2]:
        param = ParsingData('./src/interface/topo.yaml')

        fft_file = './src/interface/fft_topo/fft_' + str(i) + '.fft'
        param.set_new_value(fft_file, 'numerical_parameters', 'fft_file')
        print(param.fft_data())
        eig_freq = EigenValueProblem1D(param, 'Co', 'Py').calculate_dispersion_along_direction()
        np.savetxt('dispersion_' + str(i) + '.txt', eig_freq)

def calculate_eigen_vectors():
    for i in range(50):
        param = ParsingData('./src/interface/topo.yaml')
        fft_file = './src/interface/fft_topo/fft_' + str(i) + '.fft'
        param.set_new_value(fft_file, 'numerical_parameters', 'fft_file')

        save_file = './k_0_' + str(i)
        param.set_new_value(save_file, 'numerical_parameters', 'output_file')
        print(param.output_file('vector'))
        do_program_1D(param, 'Co', 'Py', 1e-8)


def plot_frequency(dispersion_count, color):
    for i in range(50):
        points = np.loadtxt('freq_' + str(i) + '.txt')
        points = points[dispersion_count, 2:]
        for freq_i in range(10):
            point = points[freq_i]
            plt.scatter(i/50, point/1e9, c=color)
    plt.ylim(5, 40)
    plt.xlabel('Py/Co ratio')
    plt.ylabel('Frequency (GHz)')


def plot_dispersion(ff1, ff2):
    fig, ax = plt.subplots(1, 3, figsize=(6,3))

    for a, ff in zip((ax[0], ax[1]), (ff1, ff2)):
        disp_arr = np.transpose(np.loadtxt('dispersion_' + str(ff) + '.txt'))
        for i in range(2,7):
            a.plot(np.linspace(-0.5, 0.5, len(disp_arr[i, :])),
                     disp_arr[i, :]/1e9, c='0')
            a.set_xlabel('ka/2$\pi$')
            a.set_ylim([0, 24])

            # rect = patches.Rectangle((50-ff2/2*100,-np.pi),ff2*100,2*np.pi,edgecolor='none',facecolor='C4', alpha=0.5)
            # ax[1][0].add_patch(rect)
        ax[0].set_ylabel('Frequency (GHz)')
    ax[0].set_title('ff=' + str(ff1/50))
    ax[1].set_title('ff=' + str(ff2/50))

    for k, col in zip(range(2), ('C0', 'C1')):
        for i in range(50):
            points = np.loadtxt('src/interface/disp/freq_' + str(i) + '.txt')
            points = points[k, 2:]
            for freq_i in range(5):
                point = points[freq_i]
                ax[2].scatter(i/50, point/1e9, c=col, s=0.7)
        ax[2].axvline(ff1/50, 0, 30)
        ax[2].axvline(ff2/50, 0, 30)
        ax[2].set_ylim(0, 44)
        ax[2].set_xlabel('Filling fraction')
        ax[2].set_title('k=0, k=$\pi$/a')
        # ax[2].set_ylabel('Frequency (GHz)')


    plt.savefig('dis.svg')


def plot_eigen_profile():
    ff1 = 0.3
    ff2 = 0.7
    # for i in range(15,16):
    param = ParsingData('./src/interface/topo.yaml')
    a1 = Profile1D(1, './src/interface/vector/k_0_' + str(int(ff1/2*100)) + '.vec', 'd', param)
    a2 = Profile1D(1, './src/interface/vector/k_pi_' + str(int(ff1/2*100)) + '.vec', 'd', param)
    a3 = Profile1D(1, './src/interface/vector/k_0_' + str(int(ff2/2*100)) + '.vec', 'd', param)
    a4 = Profile1D(1, './src/interface/vector/k_pi_' + str(int(ff2/2*100)) + '.vec', 'd', param)
    fig, ax = plt.subplots(2,2, figsize=(5,3))
        # for mod_num in range(5):

    x1, m1 = a1.spatial_distribution_dynamic_magnetization(5000, 1)
    print(m1.shape)
    m1 = np.roll(m1, (5000-int(ff1/2*100)*50))
    ax[0][0].plot(x1,abs(m1), c='C0')
    ax[0][0].plot(x1, np.arctan2(m1.imag, m1.real), c='C0', ls='--')
    rect = patches.Rectangle((50-ff1/2*100,-np.pi),ff1*100,2*np.pi,edgecolor='none',facecolor='C4', alpha=0.5)
    ax[0][0].add_patch(rect)
    ax[0][0].axvline(50, -np.pi, np.pi,c='C3')
    ax[0][0].set_title("k=0, ff=0.3")
    ax[0][0].set_ylabel('Amplitude / phase')


    x2, m2 = a2.spatial_distribution_dynamic_magnetization(5000, 1)
    m2 = np.roll(m2, (5000-int(ff1/2*100)*50))
    ax[1][0].plot(x2,abs(m2), c='C1')
    ax[1][0].plot(x2, np.arctan2(m2.imag, m2.real), c='C1', ls='--')
    rect = patches.Rectangle((50-ff1/2*100,-np.pi),ff1*100,2*np.pi,edgecolor='none',facecolor='C4', alpha=0.5)
    ax[1][0].add_patch(rect)
    ax[1][0].axvline(50, -np.pi, np.pi,c='C3')
    ax[1][0].set_title("k=$\pi$/a, ff=0.3")
    ax[1][0].set_ylabel('Amplitude / phase')
    ax[1][0].set_xlabel('Position (nm)')


    x3, m3 = a3.spatial_distribution_dynamic_magnetization(5000, 1)
    m3 = np.roll(m3, (5000-int(ff2/2*100)*50))
    ax[0][1].plot(x3,abs(m3), c='C0')
    ax[0][1].plot(x3, np.arctan2(m3.imag, m3.real), c='C0', ls='--')
    ax[0][1].axvline(50, -np.pi, np.pi, c='C3')
    rect = patches.Rectangle((50-ff2/2*100,-np.pi),ff2*100,2*np.pi,edgecolor='none',facecolor='C4', alpha=0.5)
    ax[0][1].add_patch(rect)
    ax[0][1].set_title("k=0, ff=0.8")

    x4, m4 = a4.spatial_distribution_dynamic_magnetization(5000, 1)
    m4 = np.roll(m4, (5000-int(ff2/2*100)*50))
    ax[1][1].plot(x4,abs(m4), c='C1')
    ax[1][1].plot(x4, np.arctan2(m4.imag, m4.real), c='C1', ls='--')
    rect = patches.Rectangle((50-ff2/2*100,-np.pi),ff2*100,2*np.pi,edgecolor='none',facecolor='C4', alpha=0.5)
    ax[1][1].add_patch(rect)
    ax[1][1].axvline(50, -np.pi, np.pi, c='C3')
    ax[1][1].set_title("k=$\pi$/a, ff=0.8")
    ax[1][1].set_xlabel('Position (nm)')
    for i in ax:
        for j in i:
            j.set_ylim(-np.pi, np.pi)
    plt.ylim(-np.pi, np.pi)
    plt.savefig('4.png', dpi=200)
    plt.clf()
    plt.close()


def plot_eigen_profile1():
    ff1 = 0.2
    ff2 = 0.5
    ff3 = 0.8
    # for i in range(15,16):
    param = ParsingData('./src/interface/topo.yaml')
    a1 = Profile1D(1, './src/interface/vector/k_pi_' + str(int(ff1/2*100)) + '.vec', 'd', param)
    a2 = Profile1D(1, './src/interface/vector/k_pi_' + str(int(ff1/2*100)) + '.vec', 'd', param)
    a3 = Profile1D(1, './src/interface/vector/k_pi_' + str(int(ff3/2*100)) + '.vec', 'd', param)
    fig, ax = plt.subplots(1,3, figsize=(6,1.5))
        # for mod_num in range(5):

    x1, m1 = a1.spatial_distribution_dynamic_magnetization(5000, 2)
    print(m1.shape)
    m1 = np.roll(m1, (5000-int(ff1/2*100)*50))
    ax[0].plot(x1,abs(m1), c='C1')
    ax[0].plot(x1, np.arctan2(m1.imag, m1.real), c='C1', ls='--')
    rect = patches.Rectangle((50-ff1/2*100,-np.pi),ff1*100,2*np.pi,edgecolor='none',facecolor='C4', alpha=0.5)
    ax[0].add_patch(rect)
    ax[0].axvline(50, -np.pi, np.pi,c='C3')
    ax[0].set_title("k=$\pi$/a, ff=0.2, ")
    ax[0].set_ylabel('Amplitude / phase')
    ax[0].set_xlabel('Position (nm)')


    x2, m2 = a2.spatial_distribution_dynamic_magnetization(5000, 2)
    m2 = np.roll(m2, (5000-int(ff2/2*100)*50))
    ax[1].plot(x2,abs(m2), c='C1')
    ax[1].plot(x2, np.arctan2(m2.imag, m2.real), c='C1', ls='--')
    rect = patches.Rectangle((50-ff2/2*100,-np.pi),ff2*100,2*np.pi,edgecolor='none',facecolor='C4', alpha=0.5)
    ax[1].add_patch(rect)
    ax[1].axvline(50, -np.pi, np.pi,c='C3')
    ax[1].set_title("k=$\pi$/a, ff=0.5")
    ax[1].set_xlabel('Position (nm)')

    x3, m3 = a3.spatial_distribution_dynamic_magnetization(5000, 2)
    m3 = np.roll(m3, (5000-int(ff3/2*100)*50))
    ax[2].plot(x3,abs(m3), c='C1')
    ax[2].plot(x3, np.arctan2(m3.imag, m3.real), c='C1', ls='--')
    ax[2].axvline(50, -np.pi, np.pi, c='C3')
    rect = patches.Rectangle((50-ff3/2*100,-np.pi),ff3*100,2*np.pi,edgecolor='none',facecolor='C4', alpha=0.5)
    ax[2].add_patch(rect)
    ax[2].set_title("k=$\pi$/a, ff=0.8")
    # ax[0][2].set_ylabel('Amplitude / phase')
    ax[2].set_xlabel('Position (nm)')


    for i in ax:
            i.set_ylim(-np.pi, np.pi)
    plt.ylim(-np.pi, np.pi)
    plt.savefig('3.png', dpi=200)
    plt.clf()
    plt.close()

if __name__ == "__main__":
    # calculate_fft(generate_sequences())
    # tab = np.loadtxt('./src/interface/periodic.fft').view(complex)
    # tab = np.transpose(np.loadtxt('./src/interface/fft_topo/fft_0.fft'))
    # tab = tab[0] + tab[1] * 1j
    # draw_elementary_cell(tab)
    # calculate_fft(generate_sequences())
    # calculate_dispersion()
    plot_dispersion(25, 30)
    # plot_frequency(0, color='C0')
    # plot_frequency(1, color='C1')
    #
    # plt.show()
    # calculate_eigen_vectors()
    # plot_eigen_profile()
    # plt.show()
