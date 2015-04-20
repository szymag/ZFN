__author__ = 'szymag'

import matplotlib.pyplot as plt


class Plot:
    def __init__(self, magnetyzacja_pod_plot):
        self.magnetyzacja_pod_plot = magnetyzacja_pod_plot

    def generowanie_danych_plot(self, k):
        list = []
        for ii in range(len(self.magnetyzacja_pod_plot[0])):
            list.append(self.magnetyzacja_pod_plot[0][ii][k])
        return list

    def plot(self):
        for ii in range(len(self.magnetyzacja_pod_plot)):
            plt.plot(self.generowanie_danych_plot(1), self.generowanie_danych_plot(2))
        plt.show()

