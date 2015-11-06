from numpy import linspace, pi, identity
from scipy import linalg

from src.MacierzDoZagadnienia import MacierzDoZagadnienia
from src.ParametryMaterialowe import ParametryMaterialowe


class ZagadnienieWlasne(ParametryMaterialowe):
    def __init__(self, rozmiar_macierzy_blok, gestosc_wektorow_q):
        ParametryMaterialowe.__init__(self, rozmiar_macierzy_blok)
        self.lista_wektorow_q = [((2 * pi * k / self.a), 0) for k in linspace(0, 1, gestosc_wektorow_q)]

    def utworz_macierz_M(self, wektor_q):
        macierz = MacierzDoZagadnienia(self.rozmiar_macierzy_blok).wypelnienie_macierzy(wektor_q)
        return macierz

    def zagadnienie_wlasne(self, wektor_q):
        m_jednostkowa = identity(2 * self.rozmiar_macierzy_blok)
        macierz_M = self.utworz_macierz_M(wektor_q)
        return linalg.eig(m_jednostkowa, macierz_M)[0]


q = ZagadnienieWlasne(103, 4)
print(q.zagadnienie_wlasne((0, 1)))
