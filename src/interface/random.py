import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from src.fft_from_image.FFT import FFT
from src.eig_problem.EigenValueProblem import EigenValueProblem1D
from src.fft_from_image.Sequences import Phason, Fibonacci, Random
from src.io.DataReader import ParsingData
from src.modes.MagnetizationProfile import Profile1D

input_parameters = ParsingData('./src/interface/Random.yaml')

repetition_seq = 5
structure_length = 377 # Fibonacci 14
inclusion_count = 140


def save_structure(inclusion_count, samples_count):
    for i in range(samples_count):
        r = Random(repetition_seq, structure_length, inclusion_count)
        seq = r.sequence()
        fft_seq = FFT().fft1d(seq)
        np.savetxt('random_' + str(inclusion_count) + '_' + str(i) + '.fft', fft_seq)
        np.savetxt('random_' + str(inclusion_count) + '_' + str(i) + '.pos', seq)


def draw_structure(inclusion_count, sample_number, ax):
    random = np.array(np.loadtxt('./random_' + str(inclusion_count) + '_' + str(sample_number) + '.pos'), dtype=int)
    for index, el in enumerate(random):
        if el == 0:
            ax.add_patch(Rectangle((index*91, 0), 91, 1, color='green', alpha=0.5, linewidth=0, edgecolor=None))
        else:
            ax.add_patch(Rectangle((index*91, 0), 91, 1, color='gray', alpha=0.5, linewidth=0, edgecolor=None))


if __name__ == "__main__":
    #save_structure(inclusion_count, 10)
    ax1 = plt.axes()
    draw_structure(inclusion_count, 1, ax1)
    ax1.set_xlim(0, 34307)
    plt.show()

