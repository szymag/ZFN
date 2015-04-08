__author__ = 'szymag'

from src.Magnetyzacja import Magnetyzacja

from src.SiecTrojkatna import SiecTrojkatna
from src.SiatkaPunktow import SiatkaPunktow
from src.Wykresy import Wykresy


siatka = SiatkaPunktow(5, 5, 1)
siec = SiecTrojkatna(2, 2, 8, 8)
w = Magnetyzacja(siatka, siec)

r = Wykresy(w.magnetyzacja_dla_sieci('okragly'))

r.wykres_pcolor()


