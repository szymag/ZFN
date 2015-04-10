__author__ = 'szymag'




class Wykresy:
    """
    Klasa odpowiadająca za prezentację wyników w formie graficznej
    """
    def __init__(self, magnetyzacja_dla_sieci):
        """
        :param magnetyzacja_dla_sieci: Tworzony jest obiekt, zawierający wartości magnetyzacji, dla danego punktu.
        Dane te tworzone są w klasie 'Magnetyzacja' w metodzie 'magnetyzacja_dla_sieci'
        """
        self.magnetyzacja_dla_sieci = magnetyzacja_dla_sieci

    def funkcja(self, k):
        """
        :param k: określa, dla którego z trzech zbiorów danych, współrzędnych x, y, magnetyzacji, tworzona jest tablica.
        :return: zwracana jest tablia zawierająca dany typ danych, w zależności od k:
         k = 1: współrzędne x
         k = 2: współrzędne y
         k = 3: magnetyzacja
        """
        table = []
        for ii in range(len(self.magnetyzacja_dla_sieci)):
            temp = []
            for jj in range(len(self.magnetyzacja_dla_sieci)):
                temp.append(self.magnetyzacja_dla_sieci[ii][jj][k])
            table.append(temp)
        return table

    def dane_do_wykresu(self):
        pass

    def wykres_pcolor(self):
        pass

    def plot(self):
        pass


