import numpy as np
import matplotlib.pyplot as plt

TM4Co = np.transpose(np.loadtxt('TM4_2809_CoPy.txt'))
TM5Co = np.transpose(np.loadtxt('TM5_2809_CoPy.txt'))
TM6Co = np.transpose(np.loadtxt('TM6_2809_CoPy.txt'))
TM4Ni = np.transpose(np.loadtxt('TM4_2809_NiFe.txt'))
TM5Ni = np.transpose(np.loadtxt('TM5_2809_NiFe.txt'))
TM6Ni = np.transpose(np.loadtxt('TM6_2809_NiFe.txt'))
TM6Niav = np.transpose(np.loadtxt('TM6_1369_NiFeav.txt'))
TM6Coav = np.transpose(np.loadtxt('TM6_841_CoPy_av.txt'))

def idos_struc_comp():
    plt.plot(TM4Co[1]/10e7, TM4Co[0])
    plt.plot(TM5Co[1]/10e7, TM5Co[0])
    plt.plot(TM6Co[1]/10e7, TM6Co[0])
    plt.xlim([65, 140])
    plt.show()

def av_comp():
    plt.plot(TM6Coav[1]/10e7, TM6Coav[0])
    #plt.plot(TM6Ni[1]/10e7, TM6Co[0])
    #plt.xlim([0, 140])
    plt.show()


av_comp()
# idos_struc_comp()