import matplotlib.pyplot as plt
import numpy as np

import sys

data_filename = 'test.txt'

"""
m = np.transpose(np.loadtxt(name))[1:]

num = np.arange(1, np.shape(m)[0] + 1)
np.savetxt(name, np.transpose([num, m]), fmt='%.0f  %.12e')
"""



def load_data(filename):
    try:
        return np.transpose(np.loadtxt(filename))
    except FileNotFoundError:
        print('No such file {}'.format(data_filename))
        sys.exit(1)


def plot_dispersion(data):
    for i in range(5):
        plt.plot(file[0], file[1 + i], color='r')
    plt.xlabel('wektor falowy q [m^-1]')
    plt.ylabel('f [Hz]')
    plt.title(name)
    #plt.ylim([0.7e10, 1.2e10])
    plt.show()

plot_dispersion()

def plot_num_freq():
    plt.plot(file[1], file[0], color='r')
    plt.xlabel('freq')
    plt.ylabel('num')
    plt.title(name)
    plt.show()


#plot_num_freq()
