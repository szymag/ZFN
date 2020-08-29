import matplotlib.pyplot as plt
import numpy as np
from math import sqrt
from matplotlib.colors import LogNorm
from matplotlib.ticker import MultipleLocator
from matplotlib.patches import Rectangle


class Plot:
    # TODO: Refactor ths class
    # TODO: Argument in this class seems to be not wise
    def __init__(self, number_of_disp_branch, x_lim=None, y_lim=None, name_of_output_file=None):
        self.name_of_file = name_of_output_file
        self.number_of_disp_branch = number_of_disp_branch
        self.x_lim = x_lim
        self.y_lim = y_lim

    def dispersion_relation(self, input_data):
        if isinstance(input_data, str):
            loaded_data = np.transpose(np.loadtxt(input_data))
        else:
            loaded_data = np.transpose(input_data)
        assert len(loaded_data.shape) == 2, 'loaded file must be two dimensional'

        for i in range(self.number_of_disp_branch):
            plt.plot(np.hypot(loaded_data[0], loaded_data[1]) / 10e6, loaded_data[2 + i] / 10e8, color='r')

        if self.y_lim is not None:
            plt.ylim(self.y_lim)

        if self.x_lim is not None:
            plt.xlim(self.x_lim)

        plt.grid()
        plt.grid()

        plt.xlabel('reciprocal vector k [mum^-1]')
        plt.ylabel('frequency [GHz]')

        self.show_or_save_plot()

    def bls(self, frequencies, weights, angle):
        plt.style.use('fivethirtyeight')
        plt.rcParams['axes.facecolor'] = '#F6FBFC'
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_aspect(aspect=0.17)
        # A = 1.1e-11
        # ms = 0.86e6
        # d = 30e-9
        # mu0 = 1.2566370614359173e-06
        # gamma = 176e9
        # H = 0.05 / mu0
        # heff = H
        # k = np.arange(1, 3e8, 50)
        # wh = gamma * mu0 * (heff + 2 * A * k ** 2 / (mu0 * ms))
        # wm = gamma * mu0 * ms
        # ax.plot(k/np.pi*200e-9/2, np.sqrt((wh + wm*(1-np.exp(-k*d))/(k*d)) * (wh)) / (2*np.pi*1e9), ls='--')
        x = np.hypot(frequencies[:, 0], frequencies[:, 1])
        y = frequencies[:, 2:]
        z = weights[:, 2:]
        for ind, i in enumerate(x):
            order = z[ind, :].argsort()
            ax.scatter(np.zeros(y.shape[1]) + i/np.pi*200e-9/2, y[ind, order]/ 1e9,
                        c=z[ind, order], cmap='BuGn')

        if self.y_lim is not None:
            ax.set_ylim(self.y_lim[0], self.y_lim[1])

        if self.x_lim is not None:
            plt.xlim(self.x_lim[0], self.x_lim[1])
        ax.set_xlabel(r'$Wave vector\ (\frac{2\pi}{a})$')
        ax.set_ylabel(r'$Frequency\ (GHz)$')
        ax.set_title(str(angle) + r'$^\circ$')
        plt.tight_layout()

        plt.savefig(str(angle) + 'deg.png', dpi=400)
        plt.show()

    def bls_for_given_frequency(self, frequencies, weights, frequency):
        plt.style.use('fivethirtyeight')
        plt.rcParams['axes.facecolor'] = '#F6FBFC'
        fig = plt.figure()
        ax = fig.add_subplot(111)
        #ax.set_aspect(aspect=0.83)
        plt.rcParams['axes.facecolor'] = '#440255'
        positions = np.argwhere(np.abs(frequencies[:, 2:] - frequency) < 0.1e9)
        x = weights[positions[:,0], 1]/np.pi*200e-9/2
        y = weights[positions[:,0], 0]/np.pi*200e-9/2
        color = weights[positions[:,0], positions[:,1] + 2]

        for i, j in [(1,1), (-1,1), (1,-1), (-1, -1)]:
            ax.scatter(i * x, j * y, c=color, s=8,
                    norm=LogNorm(vmax=weights[:, 2:].max(), vmin=weights[:, 2:].min()),
                    cmap='BuPu')

        if self.y_lim is not None:
            ax.set_ylim(self.y_lim[0], self.y_lim[1])

        if self.x_lim is not None:
            plt.xlim(self.x_lim[0], self.x_lim[1])
        ax.set_xlabel(r'$Wave 3vector\ (\frac{2\pi}{a})$')
        ax.set_ylabel(r'$Wave vector\ (\frac{2\pi}{a})$')
        ax.set_title(r'$' + str(np.around(frequency/1e9, 2)) + ' GHz$')
        plt.tight_layout()
        ax.set(aspect='equal')
        ax.set_xticks(np.arange(-1.5, 2, 0.5))
        plt.savefig(str(np.around(frequency/1e9, 2)) + 'ghz.png', dpi=400)
        plt.close()

    def contour_plot(self, input_data):
        X = input_data[:, 0].reshape(80, 80)
        Y = input_data[:, 1].reshape(80, 80)
        Z = input_data[:, 2].reshape(80, 80)
        plt.contour(X.reshape(sqrt(len(X)), sqrt(len(X))),
                    Y.reshape(sqrt(len(Y)), sqrt(len(Y))),
                    Z.reshape(sqrt(len(Z)), sqrt(len(Z))))
        plt.show()

    def idos(self, axis, input_data, color, alpha):
        axis.scatter(input_data / 1e9, np.arange(len(input_data)), color=color, alpha=alpha, s=10)

    def fmr_freq_function_of_magnetic_field(self, begin_of_name_file,
                                            start_number, end_number, scaling_factor_x_axis=1):
        f = plt.figure()
        ax = f.add_subplot(111)
        plt.tight_layout(pad=6.5, w_pad=5.5, h_pad=5.0)
        x_axis = np.arange(start_number, end_number + 1) / scaling_factor_x_axis

        data = self.load_and_join_frequency_for_diff_field(begin_of_name_file, start_number, end_number)

        for i in range(1, self.number_of_disp_branch):
            if i % 2 != 0:
                ax.generate_plot(x_axis, data[i])
            else:
                ax.generate_plot(x_axis, data[i], ls='--')

        if self.x_lim is not None:
            ax.set_xlim(self.x_lim)
        if self.y_lim is not None:
            ax.set_ylim(self.y_lim)

        ax.xaxis.grid()
        ax.yaxis.grid()

        ax.set_ylabel('frequency (GHz)', fontsize=22)
        ax.set_xlabel(r'external field $H_{0}$ (T)', fontsize=22)
        ax.set_title(r'FMR ($k=0$), $0^{\circ}$ ', fontsize=22)

        self.create_legend_for_fmr(ax)
        self.show_or_save_plot()

    def load_and_join_frequency_for_diff_field(self, begin_of_loaded_file_name, start_number, end_number):
        data = np.array(np.loadtxt(begin_of_loaded_file_name + str(start_number) + '.dat') / 10e8)
        for file_num in range(start_number + 1, end_number + 1):
            data = np.vstack((data, np.loadtxt(begin_of_loaded_file_name + str(file_num) + '.dat') / 10e8))
        return data.T

    def create_legend_for_fmr(self, axis):
        label_in_legend = []
        for i in range(1, self.number_of_disp_branch):
            label_in_legend.append('mode ' + str(i))
        return axis.legend(label_in_legend)

    def draw_structure(self, axis, sequence, phasons, stripe_width):
        seq = sequence
        # print(seq)
        for i in phasons:
            seq[(i + 1) % len(seq)] = 2
            seq[i] = -2

        stripe = lambda index, color, alpha: axis.add_patch(Rectangle((index * stripe_width, 0), stripe_width, 1,
                        color=color, alpha=alpha, linewidth=0, edgecolor=None))
        for index, el in enumerate(seq):
            if el == 0:
                stripe(index, '#f9f9f9', 1)
            elif el == 1:
                stripe(index, '#56B4E9', 1)
            elif el < 0:
                stripe(index, '#CC79A7', 1)
            elif el > 1:
                stripe(index, '#D55E00', 1)

    def show_or_save_plot(self):
        if self.name_of_file is None:
            plt.show()
        elif type(self.name_of_file) == str:
            plt.savefig(self.name_of_file + '.svg')
        else:
            plt.show()
            return 'wrong argument was put'


if __name__ == "__main__":
    Plot(5).dispersion_relation('/home/szymag/python/ZFN/src/interface/1.12e-07.txt')
