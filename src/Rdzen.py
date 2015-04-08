__author__ = 'szymag'


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

        wsp = []
        for ii in range(len(self.gx)):
            temp = []
            for jj in range(len(self.gy)):
                wspolczynnik_fouriera = self.wzory_wspolczynniki_fouriera(self.gx[ii][0], self.gy[jj][1],
                                                                    MoA, MoB, d, s, r)
                temp.append([self.gx[ii][0], self.gy[jj][1], wspolczynnik_fouriera])
            wsp.append(temp)
        return wsp

    def wylicz_wspolczynniki_fouriera(self, MoA, MoB, d, s, r):
        """ metoda zdefiniowana w klasach pochodnych """
        pass