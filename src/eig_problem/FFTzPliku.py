import numpy as np

from src.eig_problem.ParametryMaterialowe import ParametryMaterialowe


class FFTzPliku(ParametryMaterialowe):
    def __init__(self, rozmiar_macierzy_blok, filepath='fft1.txt'):
        ParametryMaterialowe.__init__(self, rozmiar_macierzy_blok)
        self.table = np.loadtxt(filepath).view(complex)
        self.shape = (len(self.table), len(self.table[0]))

    def vector_from_piksel_position(self):
        reci_vector = list(np.zeros(len(self.table) * len(self.table[0])))
        index = len(self.table[0])
        for i in range(len(self.table)):
            for j in range(len(self.table[0])):
                reci_vector[i + j * index] = (2 * np.pi * (i - int(len(self.table) / 2)) / self.a,
                                              2 * np.pi * (j - int(len(self.table[0]) / 2)) / self.a)
        assert reci_vector[len(self.table) * len(self.table[0]) - 1] != 0
        return reci_vector

    def normalized_coefficient(self):
        normalized_max_value = (self.MoCo - self.MoPy) * np.pi * self.r ** 2 / self.a ** 2 + self.MoPy
        max_value = np.amax(self.table)
        coefficient = list(np.zeros(len(self.table) * len(self.table[0])))
        index = len(self.table[0])
        for i in range(len(self.table)):
            for j in range(len(self.table[0])):
                coefficient[i + j * index] = abs(self.table[i][j]) / max_value * normalized_max_value
        assert max_value >= np.amax(abs(self.table))
        return coefficient

    def dict_vector_coeff(self):
        k = self.vector_from_piksel_position()
        v = self.normalized_coefficient()
        d = zip(k, v)
        # assert list(d.values())[200] == v[200], v[200]
        return dict(d)

    def vector_to_matrix(self):
        list_vector = self.vector_from_piksel_position()
        list_vector1 = [i for i in list_vector if abs(i[0]) <= 2 * int(self.shape[0] / 4 - 1) * np.pi / self.a
                        and abs(i[1]) <= 2 * int(self.shape[1] / 4 - 1) * np.pi / self.a]
        assert len(list_vector1) == int(np.sqrt(len(list_vector1))) ** 2, len(list_vector1)
        return list_vector1


if __name__ == "__main__":
    pass
