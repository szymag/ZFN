__author__ = 'szymag'

from src.Magnetyzacja import Magnetyzacja
from src.SiatkaPunktow import SiatkaPunktow
from src.SiecTrojkatna import SiecTrojkatna

from src.structure_from_fourier_coefficient.Plot import Plot

siatka = SiatkaPunktow(5, 5, 0.07)
siec = SiecTrojkatna(10, 10, 18, 18)
w = Magnetyzacja(siatka, siec)

r = Plot(w.magnetyzacja_dla_sieci('okragly'))

# r.wykres_pcolor()

r.plot()
