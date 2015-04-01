__author__ = 'szymag'

import numpy as np
import matplotlib.pyplot as plt


class Wykresy:
    def __init__(self, magnetyzacja_dla_sieci):
        self.magnetyzacja_dla_sieci = magnetyzacja_dla_sieci

    def funkcja(self, k):

        table = []
        for ii in range(len(self.magnetyzacja_dla_sieci)):
            temp = []
            for jj in range(len(self.magnetyzacja_dla_sieci)):
                temp.append(self.magnetyzacja_dla_sieci[ii][jj][k])
            table.append(temp)
        return table

    def dane_do_wykresu(self):
        return (np.array(self.funkcja(0)), np.array(self.funkcja(1)), np.array(self.funkcja(2)))

    def wykres_pcolor(self):
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

