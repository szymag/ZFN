import matplotlib.pyplot as plt
import numpy as np
from math import sqrt


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

    def bls(self, frequencies, weights):
        plt.rcParams['axes.facecolor'] = '#440255'
        x = np.hypot(frequencies[:, 0], frequencies[:, 1])
        y = frequencies[:, 2:]
        z = weights[:, 2:]
        for ind, i in enumerate(x):
            order = z[ind, :].argsort()
            plt.scatter(np.zeros(y.shape[1]) + i, y[ind, order],
                        c=z[ind, order], s=8)

        if self.y_lim is not None:
            plt.ylim(self.y_lim)

        if self.x_lim is not None:
            plt.xlim(self.x_lim)
        plt.show()

    def contour_plot(self, input_data, frequency):
        X = input_data[:, 0].reshape(80, 80)
        Y = input_data[:, 1].reshape(80, 80)
        Z = input_data[:, 2].reshape(80, 80)
        plt.contour(X.reshape(sqrt(len(X)), sqrt(len(X))),
                    Y.reshape(sqrt(len(Y)), sqrt(len(Y))),
                    Z.reshape(sqrt(len(Z)), sqrt(len(Z))))
        plt.show()

    def idos(self):
        pass

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
