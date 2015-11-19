from numpy import linspace, pi, savetxt
from scipy import linalg

from src.eig_problem.MacierzDoZagadnienia import MacierzDoZagadnienia
from src.eig_problem.ParametryMaterialowe import ParametryMaterialowe


class ZagadnienieWlasne(ParametryMaterialowe):
    def __init__(self, rozmiar_macierzy_blok, ilosc_wektorow_q):
        ParametryMaterialowe.__init__(self, rozmiar_macierzy_blok)
        self.lista_wektorow_q = [((2 * pi * k / self.a), 0.) for k in linspace(0.01, 0.49, ilosc_wektorow_q)] + \
                                [((2 * pi * k / self.a), 0.) for k in linspace(0.51, 0.99, ilosc_wektorow_q)]

    def utworz_macierz_M(self, wektor_q):
        assert type(wektor_q) == tuple, \
            'form of wektor_q is forbidden. wektor_q should be touple'
        assert len(wektor_q) == 2,\
            'form of wektor_q is forbidden. wektor_q should have two arguments'
        macierz = MacierzDoZagadnienia(self.rozmiar_macierzy_blok).wypelnienie_macierzy(wektor_q)
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
        czestosci_wlasne = [k * self.gamma * self.mu0H0 / (complex(1j) * 2 * pi) for k in wartosci_wlasne]
        czestosci_wlasne = [k for k in czestosci_wlasne if k.real > 0]
        # czestosci_wlasne = list(sorted(czestosci_wlasne, key=lambda x: abs(x)))
        czestosci_wlasne = [k.real for k in czestosci_wlasne]
        czestosci_wlasne = list(sorted(czestosci_wlasne))
        # for k in czestosci_wlasne:
        #    assert abs(k.real) > abs(1000 * k.imag), 'somthing is wrong: eigen frequency= ' + str(k)
        return czestosci_wlasne

    def wypisz_do_pliku(self):
        plik = []
        for k in self.lista_wektorow_q:
            tmp = [k[0]]
            tmp.extend(self.czestosci_wlasne(k))
            plik.append(tmp)
        # plik = list(zip(*plik))
        # plik.sort(key=lambda x: x[0])
        # plik  = list(zip(*plik))
        #plik = array(plik)


        savetxt('1.txt', plik)


q = ZagadnienieWlasne(729, 50)
q.wypisz_do_pliku()
# print(q.zagadnienie_wlasne((4.759988869075444527e+06, 0)))
#print(q.czestosci_wlasne((4.759988869075444527e+06, 0)))
