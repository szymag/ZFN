__author__ = 'szymag'

from src.Magnetyzacja import Magnetyzacja
from src.SiecTrojkatna import SiecTrojkatna
from src.SiatkaPunktow import SiatkaPunktow
from src.DensityPlot import DensityPlot

siatka1 = SiatkaPunktow(5, 5, 0.1)
siec1 = SiecTrojkatna(10, 10, 5, 5)
w = Magnetyzacja(siatka1, siec1)

r = DensityPlot([w.magnetyzacja_dla_sieci('okragly')])

r.wykres_pcolor()

# print(r.plot())

