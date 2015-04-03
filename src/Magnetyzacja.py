__author__ = 'szymag'

import cmath


class Magnetyzacja(object):
    def __init__(self, siatka_punktow, siec_kwadratowa):
        self.siatka_punktow = siatka_punktow
        self.siec = siec_kwadratowa

    def magnetyzacja_w_punkcie(self, parametry, x, y):
        temp = 0
        for ii in range(len(parametry)):
            for jj in range(len(parametry)):
                temp += (parametry[ii][jj][2]
                         * (cmath.exp(complex(1j) * (parametry[ii][jj][0] * x + (parametry[ii][jj][1] * y))))).real
        return temp

    def magnetyzacja(self, parametry):
        table = self.siatka_punktow.siatka()
        for ii in range(len(table)):
            for jj in range(len(table)):
                table[ii][jj][2] = \
                    self.magnetyzacja_w_punkcie(parametry, table[ii][jj][0], table[ii][jj][1])
        return table

    def magnetyzacja_dla_sieci(self, typ):
        if typ == 'kwadratowy':
            return self.magnetyzacja(self.siec.wylicz_wspolczynniki('kwadratowy'))
        elif typ == 'okragly':
            return self.magnetyzacja(self.siec.wylicz_wspolczynniki('okragly'))
        elif typ == 3:
            return None
        elif typ == 4:
            return None
        else:
            return None