__author__ = 'szymag'

import numpy as np
import matplotlib.pyplot as plt


class DensityPlot:
    """
    Klasa odpowiadająca za prezentację wyników w formie graficznej
    """

    def __init__(self, lista_magnetyzacja_dla_sieci):
        """
        :param lista_magnetyzacja_dla_sieci: Tworzony jest obiekt, zawierający wartości magnetyzacji, dla danego punktu.
        Dane te tworzone są w klasie 'Magnetyzacja' w metodzie 'magnetyzacja_dla_sieci'
        """
        self.lista_magnetyzacja_dla_sieci = lista_magnetyzacja_dla_sieci

    def generowanie_danych_density_plot(self, k):
        """
        :param k: określa, dla którego z trzech zbiorów danych, współrzędnych x, y, magnetyzacji, tworzona jest tablica.
        :return: zwracana jest tablia zawierająca dany typ danych, w zależności od k:
         k = 1: współrzędne x
         k = 2: współrzędne y
         k = 3: magnetyzacja
        """
        magnetyzacja_dla_sieci = self.lista_magnetyzacja_dla_sieci[0]
        table = []
        for ii in range(len(magnetyzacja_dla_sieci)):
            temp = []
            for jj in range(len(magnetyzacja_dla_sieci)):
                temp.append(magnetyzacja_dla_sieci[ii][jj][k])
            table.append(temp)
        return table

    def dane_do_wykresu(self):
        """
        :return: zwraca krotkę, zawierającą wszystkie tablice w typie array, wygenerowane przez metodę 'funkcja'
        """
        return np.array(self.generowanie_danych_density_plot(0)), np.array(
            self.generowanie_danych_density_plot(1)), np.array(
            self.generowanie_danych_density_plot(2))

    def wykres_pcolor(self):
        """
        :return: rysowany jest wykres na podstawie danych otrzymanych z metody 'dane_do_wykresu'
        """
        dane = self.dane_do_wykresu()
        x = dane[0]
        y = dane[1]
        z = dane[2]

        z_min, z_max = np.abs(z).min(), np.abs(z).max()

        plt.pcolor(x, y, z, cmap='RdBu', vmin=z_min, vmax=z_max)
        plt.title('magnetyzacja')
        plt.axis([x.min(), x.max(), y.min(), y.max()])
        plt.colorbar()

        return plt.show()


