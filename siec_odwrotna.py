__author__ = 'szymag'

from src.Magnetyzacja import Magnetyzacja
from src.SiecKwadratowa import SiecKwadratowa
from src.SiatkaPunktow import SiatkaPunktow


siatka = SiatkaPunktow(5, 5, 0.1)
siec = SiecKwadratowa(10, 10, 90, 5, 5)

w = Magnetyzacja(siatka, siec)
print(w.magnetyzacja_dla_sieci(1))