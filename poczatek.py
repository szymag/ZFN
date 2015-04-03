import math
from cmath import exp

import numpy as np
import scipy.special
import matplotlib.pyplot as plt


def wspolczynniki(zakres_x, zakres_y):
    """
    :type zakres_x: współrzędna x-owa wektora G w kwadratowej sieci odwrotnej
    :type zakres_y: współrzędna y-owa wektora G w kwadratowej sieci odwrotnej
    """
    nx = list(range(-zakres_x, zakres_x + 1))
    ny = list(range(-zakres_y, zakres_y + 1))
    gx = [k * (2 * math.pi / a) for k in nx]
    gy = [k * (2 * math.pi / a) for k in ny]
    wspolczynniki_wektor = []
    for ii in range(len(gx)):
        temp = []
        for jj in range(len(gy)):
            wspolczynnik_fouriera = (MoA - MoB) * d ** 2 / s * \
                                    np.sinc(np.array(gx[ii] * d / 2 / math.pi)) * np.sinc(
                np.array(gy[jj] * d / 2 / math.pi))
            temp.append([gx[ii], gy[jj], wspolczynnik_fouriera])
        wspolczynniki_wektor.append(temp)
    wspolczynniki_wektor[int(len(gx) / 2)][int(len(gy) / 2)][2] = MoA * (d ** 2 / s) + MoB * (1 - d ** 2 / s)
    return wspolczynniki_wektor


def wspolczynniki_rdzen_okrogly(zakres_x, zakres_y):
    """
    :type zakres_x: współrzędna x-owa wektora G w kwadratowej sieci odwrotnej
    :type zakres_y: współrzędna y-owa wektora G w kwadratowej sieci odwrotnej
    """
    nx = list(range(-zakres_x, zakres_x + 1))
    ny = list(range(-zakres_y, zakres_y + 1))
    gx = [k * (2 * math.pi / a) for k in nx]
    gy = [k * (2 * math.pi / a) for k in ny]
    wspolczynniki_wektor = []
    for ii in range(len(gx)):
        temp = []
        for jj in range(len(gy)):
            wspolczynnik_fouriera = 2 * (MoA - MoB) * math.pi * R ** 2 / s * \
                                    scipy.special.j1(math.sqrt(gx[ii] ** 2 + gy[jj] ** 2) * R) / (
                                        math.sqrt(gx[ii] ** 2 + gy[jj] ** 2) * R)
            temp.append([gx[ii], gy[jj], wspolczynnik_fouriera])
        wspolczynniki_wektor.append(temp)
    wspolczynniki_wektor[int(len(gx) / 2)][int(len(gy) / 2)][2] = (MoA - MoB) * math.pi * R ** 2 / s + MoB
    return wspolczynniki_wektor


def magnetyzacja_w_punkcie(wspolzedna_x, wspolzedna_y, ilosc_wektorow_x, ilosc_wektorow_y):
    """
    :type wspolzedna_x: współrzędna x-owa, dla której obliczana jest magnetyzacja
    :type wspolzedna_y: współrzędna y-owa, dla której obliczana jest magnetyzacja
    """
    temp = 0
    wspolczynniki_wektor = wspolczynniki_rdzen_okrogly(ilosc_wektorow_x, ilosc_wektorow_y)
    for ii in range(len(wspolczynniki_wektor)):
        for jj in range(len(wspolczynniki_wektor)):
            temp += (wspolczynniki_wektor[ii][jj][2] * (exp(complex(1j) * (
                (wspolczynniki_wektor[ii][jj][0] * wspolzedna_x) + (
                    wspolczynniki_wektor[ii][jj][1] * wspolzedna_y))))).real
    return temp


def generowanie_punktow(zakres, krok):
    temp = 0
    lista = []
    while zakres >= temp:
        lista.extend([temp])
        temp += krok
    return lista


def siatka(zakres_x, zakres_y, krok):
    """
    tworzenie siatki punktów na podstawie podanych zakresów i kroku
    """
    tempx = generowanie_punktow(zakres_x, krok)
    tempy = generowanie_punktow(zakres_y, krok)

    tablica = []
    for ii in tempx:
        temp = []
        for jj in tempy:
            temp.append([ii, jj, 0])
        tablica.append(temp)
    return tablica


def obliczanie_magnetyzacji(krok, ilosc_wektorow_x, ilosc_wektorow_y):
    tablica = siatka(a / 2, a / 2, krok)
    for ii in range(len(tablica)):
        for jj in range(len(tablica)):
            tablica[ii][jj][2] = magnetyzacja_w_punkcie(tablica[ii][jj][0], tablica[ii][jj][1], ilosc_wektorow_x,
                                                        ilosc_wektorow_y)
    return tablica


def drukowanie(tablica):
    """
    :type tablica: dowolna tablica składająca się z listy list
    """
    for ii in range(len(tablica)):
        print(tablica[ii])


def dane_do_wykresu(krok, ile_wektorow_x, ile_wektorow_y):
    X = np.arange(0, a / 2, krok)
    Y = np.arange(0, a / 2, krok)
    X, Y = np.meshgrid(X, Y)
    Z = []
    for ii in range(len(X)):
        temp = []
        for jj in range(len(Y)):
            temp.append(magnetyzacja_w_punkcie(X[ii][jj], Y[ii][jj], ile_wektorow_x, ile_wektorow_y))
        Z.append(temp)
    return (X, Y, np.array(Z))


def wykres_pcolor(krok, ile_wektorow_x, ile_wektorow_y):
    dane = dane_do_wykresu(krok, ile_wektorow_x, ile_wektorow_y)
    x = dane[0]
    y = dane[1]
    z = dane[2]

    z_min, z_max = np.abs(z).min(), np.abs(z).max()

    plt.pcolor(x, y, z, cmap='RdBu', vmin=z_min, vmax=z_max)
    plt.title('magnetyzacja')
    plt.axis([x.min(), x.max(), y.min(), y.max()])
    plt.colorbar()

    return plt.show()

# Definicja stałych:
# szerokość rdzenia - d (w notatkach jako r)
# szerokość komórki elementarnej
# powierzchnia; kwadrat o szerokości komórki elementarnej - s
# magnetyzacja rdzenia - MoA
# magnetyzacja wypełnienia - MoB
d = 6
a = 10
s = a ** 2
MoA = 10
MoB = 4
R = 3

wykres_pcolor(1, 2, 2)