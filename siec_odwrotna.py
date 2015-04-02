__author__ = 'szymag'

from src.Magnetyzacja import Magnetyzacja
from src.SiecKwadratowa import SiecKwadratowa
from src.SiatkaPunktow import SiatkaPunktow
from src.Wykresy import Wykresy


siatka = SiatkaPunktow(5, 5, 1)
siec = SiecKwadratowa(10, 10, 90, 10, 10)

w = Magnetyzacja(siatka, siec)
# print(w.magnetyzacja_dla_sieci(1))

r = Wykresy(w.magnetyzacja_dla_sieci(1))

r.wykres_pcolor()