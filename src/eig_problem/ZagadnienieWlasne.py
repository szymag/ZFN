from scipy import linalg
from numpy import linspace, pi, savetxt
from src.eig_problem.MacierzDoZagadnienia import MacierzDoZagadnienia
from src.eig_problem.ParametryMaterialowe import ParametryMaterialowe


class ZagadnienieWlasne(ParametryMaterialowe):
    def __init__(self, ilosc_wektorow, ilosc_wektorow_q):
        ParametryMaterialowe.__init__(self, ilosc_wektorow)
        self.lista_wektorow_q = [((2 * pi * k / self.a), 0.) for k in linspace(0.01, 0.5, ilosc_wektorow_q)]

    def utworz_macierz_M(self, wektor_q):
        assert type(wektor_q) == tuple, \
            'form of wektor_q is forbidden. wektor_q should be touple'
        assert len(wektor_q) == 2,\
            'form of wektor_q is forbidden. wektor_q should have two arguments'
        macierz = MacierzDoZagadnienia(self.ilosc_wektorow).wypelnienie_macierzy(wektor_q)
        return macierz

    def zagadnienie_wlasne(self, wektor_q):
        # m_jednostkowa = identity(2 * self.rozmiar_macierzy_blok)
        assert type(wektor_q) == tuple, \
            'form of wektor_q is forbidden. wektor_q should be touple'
        assert len(wektor_q) == 2,\
            'form of wektor_q is forbidden. wektor_q should have two arguments'
        macierz_M = self.utworz_macierz_M(wektor_q)
        return linalg.eig(macierz_M, right=False)

    def czestosci_wlasne(self, wektor_q):
        assert type(wektor_q) == tuple, \
            'form of wektor_q is forbidden. wektor_q should be touple'
        assert len(wektor_q) == 2,\
            'form of wektor_q is forbidden. wektor_q should have two arguments'
        wartosci_wlasne = self.zagadnienie_wlasne(wektor_q)
        czestosci_wlasne = [i.imag * self.gamma * self.mu0H0 / 2 / pi for i in wartosci_wlasne]
        czestosci_wlasne = [i for i in czestosci_wlasne if i > 0]
        czestosci_wlasne = list(sorted(czestosci_wlasne))
        return czestosci_wlasne

    def wypisz_do_pliku(self):
        plik = []
        for k in self.lista_wektorow_q:
            tmp = [k[0]]
            tmp.extend(self.czestosci_wlasne(k))
            plik.append(tmp)
        savetxt('1.txt', plik)
