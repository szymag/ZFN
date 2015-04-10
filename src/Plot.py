__author__ = 'szymag'

from pylab import *

from src.Wykresy import Wykresy


class Plot(Wykresy):
    def __init__(self, magnetyzacja_dla_sieci):
        Wykresy.__init__(self, magnetyzacja_dla_sieci)
        self.magnetyzacja_dla_sieci = magnetyzacja_dla_sieci

    def dane_do_wykresu(self):
        return (Wykresy.funkcja(self, 1)[0], Wykresy.funkcja(self, 2)[0])

    def plot(self):
        plot(self.dane_do_wykresu()[0], self.dane_do_wykresu()[1])
        show()