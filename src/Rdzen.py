__author__ = 'szymag'

from itertools import product

from numpy import zeros


class Rdzen:
    """
    Klasa. w której obliczane są współczynniki fourierea, dla różnych typów sieci.
    """

    def __init__(self, lista_wektorow_b1, lista_wektorow_b2):
        """
        :param lista_wektorow_b1: definiowanie listy wektorów sieci odwrtonej wzdłuż osi OX
        :param lista_wektorow_b2: definiowanie kolejnej listy wektorów sieci odwrotnej
        """
        self.gx = lista_wektorow_b1[:]
        self.gy = lista_wektorow_b2[:]

    def wzory_wspolczynniki_fouriera(self, wektor_g1, wektor_g2, MoA, MoB, d, s, r):
        """ metoda zdefiniowana w klasach pochodnych """
        pass

    def tablica_wspolczynniki_fouriera(self, MoA, MoB, d, s, r):
        tablica = zeros((len(self.gy), len(self.gx)), dtype=object)
        for ii, jj in product(range(len(self.gx)), range(len(self.gy))):
            wspolczynnik_fouriera = self.wzory_wspolczynniki_fouriera(self.gx[ii][0], self.gy[jj][1],
                                                                      MoA, MoB, d, s, r)
            tablica[ii][jj] = [self.gx[ii][0], self.gy[jj][1], wspolczynnik_fouriera]
        return tablica

    def wylicz_wspolczynniki_fouriera(self, MoA, MoB, d, s, r):
        """ metoda zdefiniowana w klasach pochodnych """
        pass