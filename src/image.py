import glob

from PIL import Image


class TablicaWartosciPikseli:
    def __init__(self):
        self.lista_plikow = list((glob.glob("*.png")))

    def wez_piksele(self):
        piksele = []
        for plik in self.lista_plikow:
            im = Image.open(plik)
            r = list(im.getdata(0))
            piksele.append(r)
        return piksele

    def stworz_tablice(self):
        pass

    def konwersja_dwa_kolory(self):
        pass


q = TablicaWartosciPikseli()

print(q.wez_piksele())
