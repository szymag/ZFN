import matplotlib.pyplot as plt
import numpy as np
file = np.transpose(np.loadtxt('225.txt'))

def plot():
    for i in range(17):
        plt.plot(file[0], file[i+1])
    plt.show()

plot()