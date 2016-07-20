import matplotlib.pyplot as plt
import numpy as np

name = 'f8000.txt'

file = np.transpose(np.loadtxt(name))

def plot_dispersion():
    for i in range(3):
        plt.plot(file[0], file[1 + i], color='r')
    plt.xlabel('wektor falowy q [m^-1]')
    plt.ylabel('f [Hz]')
    plt.title(name)
    #plt.ylim([0.7e10, 1.2e10])
    plt.show()

#plot_dispersion()

