import numpy as np
import matplotlib.pyplot as plt


class Reconstruction:
    def __init__(self, file):
        self.file = np.loadtxt(file).view(complex)
        self.size = (len(self.file[0]), len(self.file))

    def inverse_transform(self):
        return np.fft.fft2(self.file)

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
        return np.array(np.meshgrid(np.arange(0, self.size[0]), np.arange(0, self.size[1]))).T

    def tmp(self, position):
        rec_vector = self.reciprocal_vector()
        print(rec_vector.shape)

    def dft(self):
        pass

if __name__ == "__main__":
    q = ReconstructionAngleRotation('stripes.txt', 45)
    q.tmp(1)
