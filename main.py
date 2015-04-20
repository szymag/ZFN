__author__ = 'szymag'

from src.Magnetyzacja import Magnetyzacja
from src.SiecTrojkatna import SiecTrojkatna
from src.SiatkaPunktow import SiatkaPunktow
from src.Plot import Plot


siatka1 = SiatkaPunktow(5, 5, 0.01)
siec1 = SiecTrojkatna(10, 10, 30, 30)
w = Magnetyzacja(siatka1, siec1)

r = Plot([w.magnetyzacja_pod_plot('okragly')])

# r.wykres_pcolor()

print(r.plot())


