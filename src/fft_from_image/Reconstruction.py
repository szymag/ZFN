import numpy as np
import matplotlib.pyplot as plt
from math import cos, sin

class Reconstruction:
    def __init__(self, file):
        self.coefficient = np.loadtxt(file).view(complex)
        self.size = (len(self.coefficient[0]), len(self.coefficient))

    def inverse_transform(self):
        return np.fft.fftshift(np.fft.fft2(self.coefficient))

    def plot(self):
        y, x = np.mgrid[slice(0, self.size[1], 1), slice(0, self.size[0], 1)]
        z = abs(self.inverse_transform())
        plt.pcolor(x, y, z, cmap='gray')
        plt.colorbar()
        plt.show()


class ReconstructionAngleRotation(Reconstruction):
    def __init__(self, file, angle):
        Reconstruction.__init__(self, file)
        self.angle = angle

    def reciprocal_vector(self):
        co = cos(2 * np.pi * self.angle / 360)
        si = sin(2 * np.pi * self.angle / 360)
        x = 2 * np.pi * np.arange(0, self.size[0]) / 100
        y = 2 * np.pi * np.arange(0, self.size[1]) / 100
        return np.array(np.meshgrid(x * co - y * si,
                        x * si + y * co)).T

    def dft_single_point(self, position):
        rec_vector = self.reciprocal_vector()
        return np.sum(self.coefficient *
                        np.exp(1j * (rec_vector[:, :, 0] * position[0] + rec_vector[:, :, 1] * position[1])))

    def dft(self):
        x = np.linspace(0, 200, 100)
        y = np.linspace(0, 100, 100)
        values = np.zeros((len(x), len(y)), dtype=complex)
        for i in enumerate(x):
            for j in enumerate(y):
                values[i[0], j[0]] = self.dft_single_point([i[1],j[1]])
        return x, y, abs(values)

    def plot(self):
        tmp = self.dft()
        x, y = np.meshgrid(tmp[0], tmp[1])
        z = tmp[2]
        plt.pcolor(x, y, z, cmap='gray')
        plt.colorbar()
        plt.show()

if __name__ == "__main__":
    q = ReconstructionAngleRotation('stripes.txt', 90)
    q.plot()
