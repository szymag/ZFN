import glob
import os
from src.eig_problem.ParametryMaterialowe import ParametryMaterialowe
from PIL import Image
from numpy import zeros


class TablicaWartosciPikseli(ParametryMaterialowe):
    """
    Klasa konwertująca pliki graficzne '*.png' na tablice. Klasa przygotowana jest na obrazki dwukolorowe. Biały kolor
    oznacza wypełnienie, a czarny rdzeń.
    """
    def __init__(self,ilosc_wektorow='', typ_pola_wymiany='', start_path="."):
        self.lista_plikow = list((glob.glob(os.path.join(start_path, "*.png"))))
        ParametryMaterialowe.__init__(self, ilosc_wektorow, typ_pola_wymiany)

    def stworz_tablice(self, wczytany_plik):
        """
        Metoda tworzy tablice, odpowiadającą każdemu plikowi '*.png' w katalogu. Kolorowi czarnemu odpowiada wartość 1,
        a białemu 0.
        :param wczytany_plik:
        :return: Tablica wartości pikseli.
        """
        plik = Image.open(wczytany_plik)
        piksele = [abs((k / 255) - 1) for k in list(plik.getdata(0))]
        rozmiar = plik.size
        tablica_pikseli = list(zeros(rozmiar[1]))
        for i in range(0, rozmiar[1]):
            tablica_pikseli[i] = piksele[rozmiar[0] * i:(i + 1) * rozmiar[0]]

        return tablica_pikseli

    def tablica_dla_plikow(self):
        """
        Wywołanie metody 'stworz_tablice' dla wszystkich tablic przetwarzanych w klasie.
        :return: Lista tablic.
        """
        return [self.stworz_tablice(str(k)) for k in self.lista_plikow]
