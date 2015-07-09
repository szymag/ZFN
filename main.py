__author__ = 'szymag'

from src.Magnetyzacja import Magnetyzacja
from src.SiecKwadratowa import SiecKwadratowa
from src.SiatkaPunktow import SiatkaPunktow
from src.Plot import Plot

siatka1 = SiatkaPunktow(5, 5, 0.1)
siec1 = SiecKwadratowa(10, 10, 5, 5)
w = Magnetyzacja(siatka1, siec1)

siatka2 = SiatkaPunktow(5, 5, 0.1)
siec2 = SiecKwadratowa(10, 10, 10, 10)
ww = Magnetyzacja(siatka1, siec1)

r = Plot(ww.magnetyzacja_dla_sieci('kwadratowy'))


# r.wykres_pcolor()

print(r.plot())
