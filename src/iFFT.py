import matplotlib.pyplot as plt
from PIL import Image
from numpy import fft, zeros, mgrid


class FFT:
    def __init__(self, plik):
        self.plik = Image.open(plik)

    def wez_piksele(self):
        piksele = [abs((k / 255) - 1) for k in list(self.plik.getdata(0))]
        rozmiar = self.plik.size
        tablica_pikseli = list(zeros(rozmiar[1]))
        for i in range(0, rozmiar[1]):
            tablica_pikseli[i] = piksele[rozmiar[0] * i:(i + 1) * rozmiar[0]]
        return tablica_pikseli

    def fft(self):
        tablica = self.wez_piksele()
        return fft.fftshift(abs(fft.fft2(tablica, norm="ortho") / 40))

    def wykres(self):
        rozmiar = self.plik.size
        x, y = mgrid[slice(0, rozmiar[0], 1), slice(0, rozmiar[0], 1)]
        z = self.fft()
        #        z_min, z_max = -abs(z).max(), abs(z).max()
        plt.pcolor(x, y, z, cmap='RdBu')
        plt.colorbar()
        plt.show()


q = FFT("3fft.png")
print(q.fft())
q.wykres()
