import numpy as np
from src.eig_problem.ParametryMaterialowe import ParametryMaterialowe


class FFTfromFile(ParametryMaterialowe):
    """
    Klasa, która wczytuje zadany plik tekstowy (domyślnie jest to 'fft1.txt' i wyciąga z niego informację o wektorach
    sieci odrwotnej wraz z odpowiadającymi im współczynnikami Fouriera.
    """
    def __init__(self, ilosc_wektorow, coeff_number=None, filepath='fft1.txt'):
        ParametryMaterialowe.__init__(self, ilosc_wektorow)
        self.table = np.loadtxt(filepath).view(complex)
        if coeff_number is None:
            self.coeff_number = ilosc_wektorow
        else:
            assert coeff_number < ilosc_wektorow, 'matrix of element is to small to that number'
            self.coeff_number = coeff_number

    def vector_from_piksel_position(self):
        """
        Dla każdego piksela, na podstawie jego położenia, wyznaczany jest wektor sieci odwrotnej.
        :return: Lista wektorów sieci odwrotnej.
        """
        reci_vector = list(np.zeros(len(self.table) * len(self.table[0])))
        index = len(self.table[0])
        for i in range(len(self.table)):
            for j in range(len(self.table[0])):
                reci_vector[i + j * index] = (2 * np.pi * (i - int(len(self.table) / 2)) / self.a,
                                              2 * np.pi * (j - int(len(self.table[0]) / 2)) / self.b)
        assert reci_vector[len(self.table) * len(self.table[0]) - 1] != 0
        return reci_vector

    def coefficient(self):
        """
        Metoda tworząca listę współczynników. Na podstawie położenia w tablicy, określane jest położenie w liście
        :return: Lista współczynników.
        """
        coefficient = list(np.zeros(len(self.table) * len(self.table[0])))
        index = len(self.table[0])
        for i in range(len(self.table)):
            for j in range(len(self.table[0])):
                coefficient[i + j * index] = self.table[i][j]
        return coefficient

    def dict_vector_coeff(self):
        """
        Metoda ta tworzy słownik, gdzi kluczem jest wektor sieci odwrotnej, a wartością współczynnik Fouriera.
        :return: Słownik
        """
        k = self.vector_from_piksel_position()
        v = self.coefficient()
        d = zip(k, v)
        return dict(d)

    def vector_to_matrix(self):
        """
        Zwracana jest lista wektorów, służąca dalej do zagadnienia własnego. Metoda jest tak zdefiniowana, by wybierać
        wektory z kwadratu, wokół zerowego wsółczynnika Fouriera. Ilość jest tak dobrana, by pasowała do rozmiaru tablicy.
        :return: Lista wektorów, o zadanej z zewnątrz liczbie elementów.
        """
        list_vector = self.vector_from_piksel_position()
        list_vector = [i for i in list_vector if abs(i[0]) <= 2 * int(np.sqrt(self.coeff_number) / 2) * np.pi / self.a
                       and abs(i[1]) <= 2 * int(np.sqrt(self.coeff_number) / 2) * np.pi / self.b]
        assert len(list_vector) == int(np.sqrt(len(list_vector))) ** 2, len(list_vector)
        return list_vector


if __name__ == "__main__":
    pass
