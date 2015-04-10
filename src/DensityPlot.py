__author__ = 'szymag'
import numpy as np
import matplotlib.pyplot as plt

from src.Wykresy import Wykresy


class DensityPlot(Wykresy):
    def __init__(self, magnetyzacja_dla_sieci):
        Wykresy.__init__(self, magnetyzacja_dla_sieci)
        self.magnetyzacja_dla_sieci = magnetyzacja_dla_sieci

    def dane_do_wykresu(self):
        """
        :return: zwraca krotkę, zawierającą wszystkie tablice w typie array, wygenerowane przez metodę 'funkcja'
        """
        return np.array(Wykresy.funkcja(self, 0)), np.array(Wykresy.funkcja(self, 1)), np.array(
            Wykresy.funkcja(self, 2))

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