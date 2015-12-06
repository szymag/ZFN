import numpy as np

from src.eig_problem.ParametryMaterialowe import ParametryMaterialowe


class FFTfromFile(ParametryMaterialowe):
    """
    Klasa, która wczytuje zadany plik tekstowy (domyślnie jest to 'fft1.txt' i wyciąga z niego informację o wektorach
    sieci odrwotnej wraz z odpowiadającymi im współczynnikami Fouriera.
    """

    def __init__(self, ilosc_wektorow, typ_pole_wymiany, filepath='fft1.txt'):
        ParametryMaterialowe.__init__(self, ilosc_wektorow, typ_pole_wymiany)
        self.table = np.loadtxt(filepath).view(complex)
        self.coeff_number = ilosc_wektorow

    @staticmethod
    def suma_roznica_wektorow(wektor_1, wektor_2, znak):
        """
        Metoda, która w zależności od znaku oblicza sumę, bądż różnicę wektorów.
        :type wektor_1: tuple
        :type wektor_2: tuple
        :type znak: str
        :param wektor_1: Pierwszy wektor do obliczenia różnicy.
        :param wektor_2: Drugi wektor do obliczenia różnicy.
        :param znak: Określa czy obliczana ma być różnica, czy suma wektorów.
        :return: suma lub rożnica wektorów
        """
        assert len(wektor_1) == 2, \
            'form of wektor_q is forbidden. wektor_1 should have two arguments'
        assert len(wektor_2) == 2, \
            'form of wektor_q is forbidden. wektor_2 should have two arguments'
        assert type(znak) == str, \
            'znak is sign between two vector. Should be string'
        assert znak == '+' or znak == '-', \
            'only - and + are permitted'

        if znak == "-":
            return tuple([k[0] - k[1] for k in zip(wektor_1, wektor_2)])
        elif znak == "+":
            return tuple([k[0] + k[1] for k in zip(wektor_1, wektor_2)])

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

    def fourier_coefficient(self):
        """
        Metoda ta tworzy słownik, gdzi kluczem jest wektor sieci odwrotnej, a wartością współczynnik Fouriera.
        :return: Słownik
        """
        k = self.vector_from_piksel_position()
        v = [(self.MoCo - self.MoPy) * i for i in self.coefficient()]
        d = dict(zip(k, v))
        d[(0, 0)] = d[(0, 0)] + self.MoPy
        return d

    def exchange_length(self):
        """
        Metoda ta tworzy słownik, gdzi kluczem jest wektor sieci odwrotnej, a wartością współczynnik Fouriera.
        :return: Słownik
        """
        k = self.vector_from_piksel_position()
        v = [(self.lCo - self.lPy) * i for i in self.coefficient()]
        d = dict(zip(k, v))
        d[(0, 0)] = d[(0, 0)] + self.lPy
        return d


if __name__ == "__main__":
    pass
