__author__ = 'szymag'

from src.Rdzen import Rdzen


class RdzenKwadratowy(Rdzen):
    def __init__(self, lista_wektorow_b1, lista_wektorow_b2, typ_sieci):
        Rdzen.__init__(self, lista_wektorow_b1, lista_wektorow_b2)

    def wylicz_wspolczynniki(self, MoA, MoB, d, s, r, typ_sieci):
        """
        metoda zdefiniowana w klasach pochodnych
        :return: zwracana jest tablica zawierająca współczynniki fouriera dla rdzenia okrągłego.
        """
        if typ_sieci == 'kwadratowa':
            wsp = self.wspolczynniki_fouriera_tablica(MoA, MoB, d, s, r, 'kwadratowy')
            wsp[int(len(self.gx) / 2)][int(len(self.gy) / 2)][2] = \
                MoA * (d ** 2 / s) + MoB * (1 - d ** 2 / s)
            return wsp
        elif typ_sieci == 'trojkatna':
            wsp = self.wspolczynniki_fouriera_tablica(MoA, MoB, d, s, r, 'kwadratowy')
            wsp[int(len(self.gx) / 2)][int(len(self.gy) / 2)][2] = \
                MoA * (d ** 2 / s) + MoB * (1 - d ** 2 / s)
            return wsp

            wspolczynnik_fouriera = (MoA - MoB) * d ** 2 / s * \
                                    np.sinc(np.array(wektor_g1 * d / 2 / math.pi)) * \
                                    np.sinc(np.array(wektor_g2 * d / 2 / math.pi))