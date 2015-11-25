import numpy as np
from src.eig_problem.ParametryMaterialowe import ParametryMaterialowe
import matplotlib.pyplot as plt

class FFTfromFile(ParametryMaterialowe):
    def __init__(self, ilosc_wektorow, coeff_number=None, filepath='fft1.txt'):
        ParametryMaterialowe.__init__(self, ilosc_wektorow)
        self.table = np.loadtxt(filepath).view(complex)
        if coeff_number == None:
            self.coeff_number = ilosc_wektorow
        else:
            assert coeff_number < ilosc_wektorow, 'matrix of element is to small to that number'
            self.coeff_number = coeff_number

    def vector_from_piksel_position(self):
        reci_vector  = list(np.zeros(len(self.table) * len(self.table[0])))
        index = len(self.table[0])
        for i in range(len(self.table)):
            for j in range(len(self.table[0])):
                reci_vector[i + j * index] = (2 * np.pi * (i - int(len(self.table) / 2)) / self.a,
                                              2 * np.pi * (j - int(len(self.table[0]) / 2)) / self.a)
        assert reci_vector[len(self.table) * len(self.table[0]) - 1] != 0
        return reci_vector

    def normalized_coefficient(self):
        max_value = np.amax(self.table)
        coefficient = list(np.zeros(len(self.table) * len(self.table[0])))
        index = len(self.table[0])
        for i in range(len(self.table)):
            for j in range(len(self.table[0])):
                coefficient[i + j * index] = self.table[i][j]
        assert max_value >= np.amax(abs(self.table))
        return coefficient

    def dict_vector_coeff(self):
        k = self.vector_from_piksel_position()
        v = self.normalized_coefficient()
        d = zip(k, v)
        return dict(d)

    def vector_to_matrix(self):
        list_vector = self.vector_from_piksel_position()
        list_vector = [i for i in list_vector if abs(i[0]) <= 2 * int(np.sqrt(self.coeff_number) / 2 ) * np.pi / self.a
                        and abs(i[1]) <= 2 * int(np.sqrt(self.coeff_number) / 2 ) * np.pi / self.a]
        assert len(list_vector) == int(np.sqrt(len(list_vector))) ** 2, len(list_vector)
        return list_vector


    def test_coeff(self):
        table = self.table
        max_value = np.amax(self.table)
        normalized_max_value = (self.MoCo - self.MoPy) * np.pi * self.r ** 2 / self.a ** 2 + self.MoPy
        table = table * (self.MoCo - self.MoPy)
        print(np.amax(table))
        table = np.fft.ifft2(np.fft.fftshift(table)).real
        np.savetxt('ifft1.txt', table)
        x, y = np.mgrid[slice(0, 30, 1), slice(0, 30, 1)]
        z = table
        plt.pcolor(x, y, z, cmap='gray')
        plt.colorbar()
        plt.show()


if __name__ == "__main__":
    q = FFTfromFile(9)
    q.test_coeff()


