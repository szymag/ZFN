import matplotlib.pyplot as plt
import numpy as np
name = '9.6e-07nm_361_176.0GHz.txt'
file = np.transpose(np.loadtxt(name))

def plot():
    for i in range(30):
        plt.plot(file[0], file[1+i], color='r')
    plt.xlabel('wektor falowy q [m^-1]')
    plt.ylabel('f [Hz]')
    plt.title(name)

    plt.show()

plot()