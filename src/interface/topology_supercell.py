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


def generate_sequence(ff1=0.5, ff2=0.6, rec_vector_count=300):
    arr = np.zeros(rec_vector_count)
    arr[:int(rec_vector_count//2*ff1)] += 1
    arr[rec_vector_count//2:int(rec_vector_count//2*(1+ff2))] += 1
    return arr

def calculate_fft(sequence):
    tab =  np.fft.fftshift(np.fft.fft(sequence))  / len(sequence)
    tab = np.stack((tab.real, tab.imag), axis=-1)
    np.savetxt('./src/interface/fft_topo/fftsuper' + '.fft', tab)
    print('Files saved to ./fft_topo directory')


def draw_elementary_cell(fft_file):
    fft_file = fft_file.T
    fft_file = fft_file[0] + fft_file[1] * 1j
    # plt.plot(np.abs(fft_file))
    plt.plot(np.abs(np.fft.fft(fft_file)))
    plt.show()


def calculate_dispersion():
        param = ParsingData('./src/interface/topo.yaml')
        fft_file = './src/interface/fft_topo/fftsuper.fft'
        param.set_new_value(fft_file, 'numerical_parameters', 'fft_file')
        print(param.fft_data())
        eig_freq = EigenValueProblem1D(param, 'Co', 'Py').calculate_dispersion_along_direction()
        np.savetxt('dispersion_super.txt', eig_freq)

if __name__ == "__main__":
    # calculate_fft(generate_sequence())
    # draw_elementary_cell(np.loadtxt('./src/interface/fft_topo/fftsuper.fft'))
    calculate_dispersion()
