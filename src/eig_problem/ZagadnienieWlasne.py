from numpy import linspace, pi
from scipy import linalg

from src.eig_problem.MacierzDoZagadnienia import MacierzDoZagadnienia
from src.eig_problem.ParametryMaterialowe import ParametryMaterialowe


class ZagadnienieWlasne(ParametryMaterialowe):
    def __init__(self, rozmiar_macierzy_blok, gestosc_wektorow_q):
        ParametryMaterialowe.__init__(self, rozmiar_macierzy_blok)
        self.lista_wektorow_q = [((2 * pi * k / self.a), 0) for k in linspace(0, 1, gestosc_wektorow_q)]

    def utworz_macierz_M(self, wektor_q):
        macierz = MacierzDoZagadnienia(self.rozmiar_macierzy_blok).wypelnienie_macierzy(wektor_q)
        return macierz

    def zagadnienie_wlasne(self, wektor_q):
        # m_jednostkowa = identity(2 * self.rozmiar_macierzy_blok)
        macierz_M = self.utworz_macierz_M(wektor_q)
        return linalg.eig(macierz_M, right=False)

    def zagadnienie_wlasne1(self):
        # m_jednostkowa = identity(2 * self.rozmiar_macierzy_blok)
        macierz_M = self.utworz_macierz_M(self.lista_wektorow_q[20])
        return linalg.eig(macierz_M, right=False)

    def czestosciwlasne(self):
        wartosci_wlasne = self.zagadnienie_wlasne1()
        return [k * self.gamma * self.mu0H0 * complex(-1j) for k in wartosci_wlasne]


q = ZagadnienieWlasne(13, 40)

print(q.czestosciwlasne())
