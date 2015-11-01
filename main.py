__author__ = 'szymag'

from src.SiecTrojkatna import SiecTrojkatna

siec = SiecTrojkatna(10, 10, 10, 10)
print(siec.wylicz_wspolczynniki_fouriera("kwadratowy"))

'''
siatka1 = SiatkaPunktow(5, 5, 0.001)
siec1 = SiecTrojkatna(10, 10, 100, 100)
w = Magnetyzacja(siatka1, siec1)

siatka2 = SiatkaPunktow(5, 5, 0.001)
siec2 = SiecTrojkatna(10, 10, 50, 50)
ww = Magnetyzacja(siatka2, siec2)

siatka3 = SiatkaPunktow(5, 5, 0.001)
siec3 = SiecTrojkatna(10, 10, 10, 10)
www = Magnetyzacja(siatka3, siec3)

siatka4 = SiatkaPunktow(5, 5, 0.001)
siec4 = SiecTrojkatna(10, 10, 5, 5)
wwww = Magnetyzacja(siatka4, siec4)

r = Plot([w.magnetyzacja_pod_plot('kwadratowy'),
          ww.magnetyzacja_pod_plot('kwadratowy'),
          www.magnetyzacja_pod_plot('kwadratowy'),
          wwww.magnetyzacja_pod_plot('kwadratowy')])
'''

# r.wykres_pcolor()

# print(r.plot())
