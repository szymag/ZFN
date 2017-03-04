import matplotlib.pyplot as plt
import numpy as np

import sys

data_filename = 'test.txt'

"""
m = np.transpose(np.loadtxt(name))[1:]

num = np.arange(1, np.shape(m)[0] + 1)
np.savetxt(name, np.transpose([num, m]), fmt='%.0f  %.12e')
"""

def usage():
    print("Usage:\n\t{} {}".format(sys.argv[0], 'd/n'))


def load_data(filename):
    try:
        return np.transpose(np.loadtxt(filename))
    except FileNotFoundError:
        print('No such file {}'.format(data_filename))
        sys.exit(1)


def plot_dispersion(data):
    for i in range(5):
        plt.plot(data[0], data[1 + i], color='r')
    plt.xlabel('wektor falowy q [m^-1]')
    plt.ylabel('f [Hz]')
    plt.title(data)
    #plt.ylim([0.7e10, 1.2e10])
    plt.show()


def plot_num_freq(data):
    plt.plot(data[1], data[0], color='r')
    plt.xlabel('freq')
    plt.ylabel('num')
    plt.title(data_filename)
    plt.show()


if __name__ == '__main__':
    data = load_data(data_filename)
    if len(sys.argv) == 2:
        if sys.argv[1] == 'd':
            plot_dispersion(data)
        elif sys.argv[1] == 'n':
            plot_num_freq(data)
        else:
            usage()
            sys.exit(1)
    else:
        usage()
        plot_dispersion(data)
