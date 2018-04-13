import numpy as np
import matplotlib.pyplot as plt

file1 = np.transpose(np.loadtxt('1tm5.txt'))
file2 = np.transpose(np.loadtxt('1tm7.txt'))
file3 = np.transpose(np.loadtxt('1tm9.txt'))
num = np.arange(1, np.shape(file1)[0] + 1)
def idos_struc_comp():
    plt.plot(file1, num)
    plt.plot(file2, num)
    plt.plot(file3, num)
    plt.show()


#av_comp()
idos_struc_comp()