import numpy as np
import matplotlib.pyplot as plt


class Reconstruction:

    def __init__(self, file):
        self.file = np.loadtxt(file).view(complex)
        self.size = (len(self.file[0]), len(self.file))

    def inverse_transform(self):
        return np.fft.fft2(self.file)

    def plot(self):
        x, y = np.mgrid[slice(0, self.size[1], 1), slice(0, self.size[0], 1)]
        z = abs(self.inverse_transform())
        print(z.shape)
        plt.pcolor(x, y, z, cmap='gray')
        plt.colorbar()
        plt.show()


if __name__ == "__main__":
    q = Reconstruction('F7.txt')
    q.plot()