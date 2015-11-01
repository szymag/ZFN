__author__ = 'szymag'

import matplotlib.pyplot as plt


class Plot:
    def __init__(self, magnetyzacja_pod_plot):
        self.magnetyzacja_pod_plot = magnetyzacja_pod_plot

    def generowanie_danych_plot(self, k, index):
        list = []
        for ii in range(len(self.magnetyzacja_pod_plot[index])):
            list.append(self.magnetyzacja_pod_plot[index][ii][k])
        return list

    def plot(self):
        for ii in range(len(self.magnetyzacja_pod_plot)):
            plt.plot(self.generowanie_danych_plot(1, ii), self.generowanie_danych_plot(2, ii))
        plt.show()
