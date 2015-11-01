__author__ = 'szymag'
"""
from src import Siec
import numpy as np
from src import SiecKwadratowa


class IFFT:

    def __init__(self, siec):
        self.siec = siec

    def fourier_coefficients(self, siec):
        print(siec.wylicz_wspolczynniki_fouriera())

    def inverse_fft(self):
        pass
"""

import numpy as np
import matplotlib.pyplot as plt

Fs = 150  # sampling rate
Ts = 1.0 / Fs  # sampling interval
t = np.arange(0, 1, Ts)  # time vector
ff = 1  # frequency of the signal
y = np.sin(2 * np.pi * ff * t)

plt.subplot(2, 1, 1)
plt.plot(t, y, 'k-')
plt.xlabel('time')
plt.ylabel('amplitude')

plt.subplot(2, 1, 2)
n = len(y)  # length of the signal
k = np.arange(n)
T = n / Fs
frq = k / T  # two sides frequency range
freq = frq[range(int(n / 2))]  # one side frequency range

Y = np.fft.fft(y) / n  # fft computing and normalization
print(len(Y))
Y = Y[range(int(n / 2))]
print(len(Y))

plt.plot(freq, abs(Y), 'r-')
plt.xlabel('freq (Hz)')
plt.ylabel('|Y(freq)|')

plt.show()
