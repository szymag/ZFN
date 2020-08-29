import numpy as np
from scipy.signal import find_peaks
from src.fft_from_image.FFT import FFT
import matplotlib.pyplot as plt
import matplotlib
from src.io.DataReader import ParsingData
from src.modes.MagnetizationProfile import Profile1D

dx_element = 20
w = 1/5
elements_count = 377

dx_count = dx_element * elements_count
x = np.arange(dx_count)

input_parameters = ParsingData('./src/interface/Rychly.yaml')
loaded_modes = '/media/szymag/Dane/ZFN/JEMS2019/fib/F_5_0.vec'
idos = np.loadtxt('/media/szymag/Dane/ZFN/JEMS2019/fib/F_5_0.dys')
defects = np.loadtxt('/media/szymag/Dane/ZFN/JEMS2019/fib/F_5_0.pos')
defects = defects * dx_element

modes = Profile1D(0, loaded_modes, None, input_parameters)




def count_distance_from_defect(defect_pos):
    a =  (np.abs(x - defect_pos))**w
    b = (dx_count - np.abs(x - defect_pos))**w
    c = np.stack((a, b))
    return np.min(c, axis=0)
    # return a/np.max(a)

def calculate_d():
    d = np.ones(dx_count)
    for i in defects:
        d = d*count_distance_from_defect(i)
    return d / np.sum(d)


def calculate_localization():
    el = 400
    lam = np.zeros(el)
    d = calculate_d()
    for idx in range(el):
        mode_1 = np.log(np.abs(modes.spatial_distribution_dynamic_magnetization(dx_count, idx)[1]))
        mode_1 = mode_1 / np.max(mode_1)
        lam[idx] = np.sum(mode_1*d)

    return lam, (1 - lam)/(1+lam)

lam, m = calculate_localization()

plt.scatter(idos, lam)
plt.scatter(idos, m)
plt.show()
