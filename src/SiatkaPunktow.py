__author__ = 'szymag'


class SiatkaPunktow(object):
    def __init__(self, zakres_x, zakres_y, krok):

        """
        :param zakres_x: opisuje dla jakiego zakresu na osi x zoastanie obliczona magetyzacja. Przy czy zakres wynosi od
        -zakres do zakres
        :param zakres_y: opisuje dla jakiego zakresu na osi y zoastanie obliczona magetyzacja. Przy czy zakres wynosi od
        -zakres do zakres
        :param krok: Zmienna dyskretyzująca przedział, dzieląca <-zakres, zakres> na 2*zakres+1 przedziałów. Dla każdego
        zostanie policzona magnetyzacja.
        """
        self.zakres_x = zakres_x
        self.zakres_y = zakres_y
        self.krok = krok

    def generowanie_punktow(self, liczba):
        temp = 0
        lista = []
        while liczba >= temp:
            lista.extend([temp])
            temp += self.krok
        return lista

    def siatka(self):

        tempx = self.generowanie_punktow(self.zakres_x)
        tempy = self.generowanie_punktow(self.zakres_y)

        tablica = []
        for ii in tempx:
            temp = []
            for jj in tempy:
                temp.append([ii, jj, 0])
            tablica.append(temp)
        return tablica
