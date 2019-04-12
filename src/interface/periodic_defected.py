import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from src.fft_from_image.FFT import FFT
from src.eig_problem.EigenValueProblem import EigenValueProblem1D
from src.fft_from_image.Sequences import Phason, Fibonacci, Random, Periodic
from src.io.DataReader import ParsingData
from src.modes.MagnetizationProfile import Profile1D

input_parameters = ParsingData('./src/interface/Periodic.yaml')

repetition_seq = 10
structure_length = 377

phasons = 0.02

def save_structure(phasons_percentage, samples_count):
    for i in range(samples_count):
        f = Phason('P', repetition_seq, structure_length, phasons_percentage)
        seq, phason_positions = f.sequence_shuffling(f.seq)

        seq = f.sequence(seq)
        fft_seq = FFT().fft1d(seq)
        np.savetxt('periodic_' + str(phasons_percentage) + '_' + str(i) + '.fft', fft_seq)
        np.savetxt('periodic_' + str(phasons_percentage) + '_' + str(i) + '.pos', phason_positions)


def draw_structure(phasons_percentage, sample_number, ax):
    phasons = np.array(np.loadtxt('./periodic_' + str(phasons_percentage) + '_' + str(sample_number) + '.pos'), dtype=int) + 1
    fib = Periodic(repetition_seq, structure_length)
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
def save_eig_vector(phasons_percentage, sample_number):
    np.savetxt('periodic_' + str(phasons_percentage) + '_' + str(sample_number) + '.vec',
               calculate_eig_vector(phasons_percentage, sample_number).view(float))

def calculate_eig_vector(phasons_percentage, sample_number):
    fft_file = './periodic_' + str(phasons_percentage) + '_' + str(sample_number) + '.fft'
    print(fft_file)
    input_parameters.set_new_value(fft_file, 'numerical_parameters', 'fft_file')
    return EigenValueProblem1D(input_parameters, 'Py', 'Co').calculate_eigen_vectors(bloch_vector=[1e4, 0])


def mode(mode_number, phasons_percentage, sample_number, ax):
    vec_name = './periodic_' + str(phasons_percentage) + '_' + str(sample_number)+'.vec'
    print(vec_name)
    mod = Profile1D(mode_number, vec_name,
                    None, input_parameters)
    return mod.generate_plot(ax, 0)


if __name__ == "__main__":
    save_structure(phasons, 10)
    # ax1 = plt.axes()
    # draw_structure(phasons, 1, ax1)
    # ax1.set_xlim(0, 34307)
    # plt.show()
    save_eig_vector(phasons, 1)
    for i in range(75, 352):
        ax1 = plt.axes()
        draw_structure(phasons, 1, ax1)
        mode(i, phasons, 1, ax1)
        plt.show()
        plt.clf()
        plt.close()
