import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
from math import radians, sin
from src.eig_problem.ParametryMaterialowe import ParametryMaterialowe
from src.eig_problem.WektorySieciOdwrotnej import WektorySieciOdwrotnej
import matplotlib.gridspec as gridspec
import matplotlib.image as mpimg
from mpl_toolkits.axes_grid1 import make_axes_locatable

from matplotlib import rcParams

rcParams.update({'figure.autolayout': True})

plt.rc('text', usetex=True)
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Helvetica']

class MidpointNormalize(colors.Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))

from multiprocessing import Pool

from src.eig_problem.InputParameter import InputParameter
from src.eig_problem.ReciprocalVector import ReciprocalVector
from src.io.DataReader import ParsingData
from src.eig_problem.LoadFFT import LoadFFT2D


class Profile2D:
    def __init__(self, grid, load_data):
        # TODO: How to pass lattice constant? InputParameter is depreciated
        self.fourier_coefficients = np.loadtxt(load_data).view(complex)
        self.lattice_const_x = InputParameter.a
        self.lattice_const_y = InputParameter.b
        self.grid = grid

    def generate_plot(self, mode_number):
        lista_x, lista_y, lista_wartosci = self.spatial_distribution_dynamic_magnetization(mode_number)
        x, y = np.meshgrid(lista_x, lista_y)
        plt.pcolor(x, y, np.array(lista_wartosci))
        plt.colorbar()
        plt.show()

    def spatial_distribution_dynamic_magnetization(self, mode_number):
        mode = self.fourier_coefficients[mode_number, :]
        mode = mode[:len(mode) // 2]
        return self.reconstruct_whole_structure(mode)

    def reconstruct_whole_structure(self, input_coefficient):
        x = np.linspace(-self.lattice_const_x, self.lattice_const_x, self.grid)
        y = np.linspace(-self.lattice_const_x, self.lattice_const_y, self.grid)
        m = np.zeros(self.grid * self.grid, dtype=complex)
        for i, j in enumerate(np.dstack(np.meshgrid(x, y)).reshape(-1, 2)):
            m[i] = self.inverse_discrete_fourier_transform(input_coefficient, j)
        return x, y, m.reshape((self.grid, self.grid))

    def inverse_discrete_fourier_transform(self, data, position):
        reciprocal_vectors = 2 * np.pi * ReciprocalVector(max(data.shape)).lista_wektorow2d('min') /\
                             np.array((self.lattice_const_x, self.lattice_const_y))
        return np.sum(data * np.prod(np.exp(1j * reciprocal_vectors * position), axis=1))

    def fmr(self, mode_number):
        mode = self.spatial_distribution_dynamic_magnetization(mode_number)[2]
        return np.sum(np.abs(np.trapz(mode, axis=1)) ** 2 / np.trapz(np.abs(mode) ** 2, axis=1))

    def fmr_intensity_order(self):
        fmr = Pool().map(self.fmr, range(10))
        return np.array(fmr)


class Profile1D:
    def __init__(self, mode_number, load_data, name_of_file, **kwargs):
        self.eig_vectors = np.loadtxt(load_data).view(complex)
        self.lattice_const = InputParameter.a
        self.mode_number = mode_number - 1
        self.name_of_file = name_of_file

        if 'angle' in kwargs:
            self.angle = kwargs['angle']
        if 'field' in kwargs:
            self.field = kwargs['field']

    def generate_plot(self, ax, color_index, dummy=False):
        colors = ['C0', 'C3', 'C2']
        x, magnetization = self.spatial_distribution_dynamic_magnetization(500, self.mode_number)
        phase = np.abs(np.arctan2(magnetization.imag, magnetization.real))
        parameter = np.array([0 if i < 0.5 * np.pi else 1 for i in phase])
        if dummy:
            magnetization1 = np.ma.masked_where(parameter == 0, abs(magnetization) ** 2)
            magnetization2 = np.ma.masked_where(parameter != 0, abs(magnetization) ** 2)
        else:
            magnetization1 = np.ma.masked_where(parameter != 0, abs(magnetization) ** 2)
            magnetization2 = np.ma.masked_where(parameter == 0, abs(magnetization) ** 2)
        ax.plot(x, magnetization1, colors[color_index] + '-',
                x, magnetization2, colors[color_index] + '--')

        #ax.set_ylim(-0.05, 3)
        plt.setp(ax, xticks=[-1100, 0, 1100], xticklabels=[r'$-\Lambda$', 0, r'$\Lambda$'],
                 yticks=[0])
        ax.set_xlabel('Position', fontsize=12)
        ax.set_ylabel(r'$\left|\textup{m}_{\textup{z}}\right|$', fontsize=12)
        # self.output_plot()

    def output_plot(self):
        if self.name_of_file is None:
            plt.show()
        elif type(self.name_of_file) == str:
            plt.savefig(self.name_of_file + '.svg')
            plt.clf()
            plt.close()
        else:
            plt.show()
            return 'wrong argument was put'

    def save_to_file(self):
        to_file = self.spatial_distribution_dynamic_magnetization(500, self.mode_number)[1].real, \
                  self.spatial_distribution_dynamic_magnetization(500, self.mode_number)[1].imag

        np.savetxt(self.name_of_file + '.txt', np.transpose(to_file))

    def spatial_distribution_dynamic_magnetization(self, grid, mode_number):
        mode = self.eig_vectors[mode_number, 0:self.eig_vectors.shape[1] // 2]
        x = np.linspace(-self.lattice_const, self.lattice_const, grid)
        tmp = np.zeros(grid, dtype=complex)
        for ind in enumerate(x):
            tmp[ind[0]] = self.inverse_discrete_fourier_transform(mode, ind[1])
        return x * 10 ** 9, tmp

    def elementary_cell_reconstruction(self, grid):
        coefficient = np.transpose(np.loadtxt('c_coef_100.txt').view(complex))
        x = np.linspace(-self.lattice_const, self.lattice_const, grid)
        tmp = np.zeros(grid)
        for ind in enumerate(x):
            tmp[ind[0]] = abs(self.inverse_discrete_fourier_transform(coefficient, ind[1]))
        return x * 10 ** 9, tmp / (10 / 7) + 0.3

    def inverse_discrete_fourier_transform(self, data, vector_position):
        reciprocal_vectors = np.array(
            2 * np.pi * WektorySieciOdwrotnej(max(data.shape)).lista_wektorow1d('min') / self.lattice_const)
        reciprocal_vectors = np.array(2 * np.pi * ReciprocalVector(max(data.shape)).lista_wektorow1d('min')
                                      / self.lattice_const)

        return np.sum(data * np.exp(1j * reciprocal_vectors * vector_position))

    def fmr_intensity(self, mode):
        return np.abs(np.trapz(mode)) ** 2 / np.trapz(np.abs(mode) ** 2)

    def fmr_intensity_order(self):
        fmr = [self.fmr_intensity(self.spatial_distribution_dynamic_magnetization(500, i)[1]) for i in range(50)]
        # np.savetxt(self.name_of_file, np.array(fmr))
        return np.array(fmr)

    def acoustic_cross_section(self, grid, data=None):
        tmp = self.spatial_distribution_dynamic_magnetization(grid, data)
        x = np.linspace(-self.lattice_const, self.lattice_const, grid)
        return sin(radians(2 * self.angle)) * np.trapz(tmp[1] * np.cos(x * 2 * np.pi / self.lattice_const))


if __name__ == "__main__":
    def modes(modulation, start, stop, step, modes_count, ax, only_max_fmr=True):
        if only_max_fmr:
            try:
                fmr_map = np.loadtxt('fmr' + str(modulation) + '.dat')
            except:
                fmr(modulation, start, stop, step)
                fmr_map = np.loadtxt('fmr' + str(modulation) + '.dat')

            for i in np.arange(start, stop, step):
                modes = np.argsort(fmr_map[:, i])[-modes_count:] + 1
                for mode in enumerate(modes[::-1]):
                    if (start == 10 and mode[0] == 1) or (start == 17 and mode[0] == 0):
                        Profile1D(mode[1], '/home/szymag/python/ZFN/src/eig_problem/Ni/' + str(int(i)) + '.dat',
                                  'deg' + str(i) + '_mode' + str(mode[1]), angle=i).generate_plot(ax, 2)
                    elif start == 17 and mode[0] == 1:
                        Profile1D(mode[1], '/home/szymag/python/ZFN/src/eig_problem/Ni/' + str(int(i)) + '.dat',
                                  'deg' + str(i) + '_mode' + str(mode[1]), angle=i).generate_plot(ax, 0, True)
                    elif start == 10 and mode[0] == 0:
                        Profile1D(mode[1], '/home/szymag/python/ZFN/src/eig_problem/Ni/' + str(int(i)) + '.dat',
                                  'deg' + str(i) + '_mode' + str(mode[1]), angle=i).generate_plot(ax, 0, True)
                    else:
                        Profile1D(mode[1], '/home/szymag/python/ZFN/src/eig_problem/Ni/' + str(int(i)) + '.dat',
                                  'deg' + str(i) + '_mode' + str(mode[1]), angle=i).generate_plot(ax, mode[0])
        else:
            for i in enumerate(np.arange(start, stop, step)):
                for j in range(5, 8):
                    Profile1D(j, '/home/szymag/python/ZFN/src/eig_problem/' + str(int(i[1])) + '.dat',
                              'deg' + str(i[1]) + '_mode' + str(j), angle=i[1]).generate_plot(ax, i[0])


    def fmr(modulation, start, stop, step):
        fmr = np.zeros((50, int((stop - start) / step)))
        for i in enumerate(np.arange(start, stop, step)):
            fmr[:, i[0]] = Profile1D(1, '/home/szymag/python/ZFN/src/eig_problem/Ni/' + str(i[1]) + '.dat', None,
                                     angle=i[1]).fmr_intensity_order()
        np.savetxt('fmr' + str(modulation) + '.dat', fmr)


    def fmr_1(modulation, field, start, stop, step):
        return Profile1D(1, '/home/szymag/python/ZFN/src/eig_problem/' + str(field)
                         + '_' + str(modulation) + '.dat', None).fmr_intensity_order()


    def fmr_intensity_map(modulation, start, stop, step, ax):
        mateusz_plot(ax)
        if type(step) is float:
            y = np.loadtxt(
                '/home/szymag/python/ZFN/src/eig_problem/densefreq_vs_angle_' + str(modulation) + '.dat') / 1e9
        else:
            y = np.loadtxt('/home/szymag/python/ZFN/src/eig_problem/Ni/freq_vs_angle_' + str(modulation) + '.dat') / 1e9

        if int((stop - start) / step) != y.shape[1]:
            raise ValueError('Current settings is different from ZagadnienieWlasne')
        x = np.arange(start, stop, step)
        try:
            if type(step) is int:
                z = np.loadtxt('fmr' + str(modulation) + '.dat')
            else:
                z = np.loadtxt('densefmr' + str(modulation) + '.dat')
            if z.shape[1] != y.shape[1]:
                fmr(modulation, start, stop, step)
                z = np.loadtxt('fmr' + str(modulation) + '.dat')
        except:
            fmr(modulation, start, stop, step)
            z = np.loadtxt('fmr' + str(modulation) + '.dat')
        max = np.max(z)
        min = np.min(z)
        drawing_mode_count = 40
        for i in enumerate(np.arange(start, stop, step)):
            order = z[0:drawing_mode_count, i[0]].argsort()
            tmp = ax.scatter(np.zeros(drawing_mode_count) + i[1],
                             y[0:drawing_mode_count, i[0]][order],
                             c=z[0:drawing_mode_count, i[0]][order],
                             s=10, edgecolors='', vmin=min, vmax=max,
                             cmap=plt.cm.Wistia, alpha=0.96, norm=MidpointNormalize(midpoint=30))
        ax.set_ylim([4, 8])
        ax.axes.get_xaxis().set_ticklabels([])
        ax.set_ylabel('Frequency (GHz)', fontsize=12)
        # part responsible for color bar
        divider = make_axes_locatable(ax)
        cax2 = divider.append_axes("right", size="3%", pad=0.0)
        a = plt.colorbar(tmp, ax=ax, ticks=[], cax=cax2)
        a.set_label('Detection intensity (a.u.) (a) (b) (c)')



    def cross_section(modulation):
        fmr_map = np.loadtxt('fmr' + str(modulation) + '.dat')
        mode_freq = np.loadtxt('/home/szymag/python/ZFN/src/eig_problem/Ni/freq_vs_angle_' + str(modulation) + '.dat')
        assert fmr_map.shape == mode_freq.shape
        points = np.zeros(fmr_map.shape[1], dtype=complex)
        for i in enumerate(np.arange(fmr_map.shape[1])):
            mode = np.argsort(fmr_map[:, i[1]])[-1]
            points[i[0]] += Profile1D(mode, '/home/szymag/python/ZFN/src/eig_problem/Ni/' + str(i[1]) + '.dat',
                                      None, angle=i[1]).acoustic_cross_section(500, mode) * mode_weight(
                mode_freq[mode, i[1]])
        fig = plt.figure()
        ax = fig.add_subplot(211)
        ax.plot(np.arange(0, 91), np.abs(points))
        np.savetxt('cross_section_Ni.txt', np.transpose(np.vstack((np.arange(0, 91), np.abs(points) / np.max(np.abs(points))))))
        ax.set_xlabel('Magnetic Field Angle (Deg)', fontsize=18)
        ax.set_ylabel('Arbitrary Unit (a.u.)', fontsize=18)
        ax.set_yticks([], [])
        #plt.savefig('cross.section.png')
        plt.show()


    def mode_weight(position):
        peak_position = 4.8e9
        fwhm = .2e9
        return 8e17 / ((peak_position - position) ** 2 + (fwhm / 2) ** 2)


    def mateusz_plot(ax):
        mateusz = np.loadtxt('test.txt', delimiter=',')
        mateusz = 2 * mateusz + 1.5 * np.roll(mateusz, 2, axis=0) + np.roll(mateusz, 4, axis=0) + \
                  np.roll(mateusz, 6, axis=0) + 1.5 * np.roll(mateusz, -2, axis=0) + np.roll(mateusz, -4, axis=0)
        mateusz = mateusz / np.max(mateusz)
        ax.imshow(mateusz, cmap=plt.cm.binary, origin='lower', aspect='auto',
                        extent=[0, 90, 0, 11111111111.1 / 1e9], interpolation='lanczos', vmin=0, vmax=0.5)
        plt.ylim([4.4, 10])


    def get_y_coordinates(angles, modulation):
        fmr_map = np.loadtxt('fmr' + str(modulation) + '.dat')
        mode_freq = np.loadtxt('/home/szymag/python/ZFN/src/eig_problem/Ni/freq_vs_angle_' + str(modulation) + '.dat')
        return [mode_freq[np.argsort(fmr_map[:, i])[-2:], i] / 1e9 for i in angles]


    def make_final_plot():
        angles = [10, 17, 20, 30, 40, 60]
        y_coordinates = get_y_coordinates(angles, 40)
        labels_1 = ['(VII)', '(VIII)', '(IX)', '(X)', '(XI)', '(XII)']
        labels_2 = ['(I)', '(II)', '(III)', '(IV)', '(V)', '(VI)']
        fig = plt.figure(figsize=(6, 7))
        fig.tight_layout()
        img = mpimg.imread('struct.png')

        outer = gridspec.GridSpec(4, 1, height_ratios=[0.3, 2, 5, 1])

        gs0 = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=outer[0],
                                               hspace=0.0,
                                               wspace=0.0)
        ax00 = plt.subplot(gs0[:,0])
        ax01 = plt.subplot(gs0[:,1])
        ax00.imshow(img)
        ax01.imshow(img)
        ax00.axes.get_xaxis().set_visible(False)
        ax00.axes.get_yaxis().set_visible(False)
        ax01.axes.get_xaxis().set_visible(False)
        ax01.axes.get_yaxis().set_visible(False)
        gs1 = gridspec.GridSpecFromSubplotSpec(3, 2, subplot_spec=outer[1],
                                               hspace=0.0,
                                               wspace=0.0)
        ax1 = plt.subplot(gs1[0, 0])
        ax2 = plt.subplot(gs1[0, 1])
        ax3 = plt.subplot(gs1[1, 0])
        ax4 = plt.subplot(gs1[1, 1])
        ax5 = plt.subplot(gs1[2, 0])
        ax6 = plt.subplot(gs1[2, 1])
        for ind, ax in enumerate([ax1, ax3, ax5, ax2, ax4, ax6]):
            modes(40, angles[ind], angles[ind] + 1, 1, 2, ax)
            ax.text(0.05, 0.7, labels_1[ind], ha="center",
                    transform=ax.transAxes)
            ax.text(0.05, 0.2, labels_2[ind], ha="center",
                    transform=ax.transAxes)

        ax1.axes.get_xaxis().set_visible(False)
        ax2.axes.get_xaxis().set_visible(False)
        ax2.axes.get_yaxis().set_visible(False)
        ax3.axes.get_xaxis().set_visible(False)
        ax4.axes.get_xaxis().set_visible(False)
        ax4.axes.get_yaxis().set_visible(False)
        ax6.axes.get_yaxis().set_visible(False)

        gs2 = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=outer[2], hspace=0.0, wspace=0.0)
        ax7 = plt.subplot(gs2[0, :])
        fmr_intensity_map(40, 0, 91, 1, ax7)
        for i, j, k in zip(angles, labels_1, y_coordinates):
            if j == '(a)':
                ax7.text(i, 1.03*k[1], j, ha="center")
            else:
                ax7.text(i, 0.92*k[1], j, ha="center")
        for i, j, k in zip(angles, labels_2, y_coordinates):
            if j == '(g)':
                ax7.text(i, 0.92*k[0], j, ha="center")
            else:
                ax7.text(i, 1.03*k[0], j, ha="center")

        gs3 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=outer[3], hspace=0.0, wspace=0.0)
        ax8 = plt.subplot(gs3[:, :])
        ax8.plot(*np.transpose(np.loadtxt('cross_section_Ni_corrected.txt')))
        ax8.set_yticks([0])
        ax8.set_ylabel('Prec. Ampl. (a.u.)', fontsize=12)
        ax8.set_xlabel('Magnet Angle (deg)', fontsize=12)

        plt.savefig('final.svg')


    def show_modes_grid():
        fields = [0.03, 0.035, 0.04, 0.045, 0.05, 0.055, 0.06, 0.065]
        modulations = [793, 837, 877, 894, 903, 909]

        grid = gridspec.GridSpec(6, 8, hspace=0.0, wspace=0.0)
        f = plt.figure()
        f.suptitle('Angle=60 Deg')
        for ind_1, modulation in enumerate(modulations):
            frequencies = np.around(
                np.loadtxt('/home/szymag/python/ZFN/src/eig_problem/freq_vs_angle_' + str(modulation) + '.dat') / 1e9,
                decimals=2)
            for ind_2, field in enumerate(fields):
                mode_number = np.argsort(fmr_1(modulation, field, 0, 1, 1))[-1] + 1
                Profile1D(mode_number, '/home/szymag/python/ZFN/src/eig_problem/'
                          + str(field) + '_' + str(modulation) + '.dat', 'dummy',
                          field=field, angle=60).generate_plot(plt.subplot(grid[ind_1, ind_2]), 0)

                plt.subplot(grid[ind_1, ind_2]).text(0.2, 0.1, r'$H_{0}=' + str(field) + '$', ha="center",
                                                     transform=plt.subplot(grid[ind_1, ind_2]).transAxes)
                plt.subplot(grid[ind_1, ind_2]).text(0.2, 0.8, r'$M_{min}=0.' + str(modulation) + '$', ha="center",
                                                     transform=plt.subplot(grid[ind_1, ind_2]).transAxes)
                plt.subplot(grid[ind_1, ind_2]).text(0.5, 0.5, r'$freq=' +
                                                     str(frequencies[mode_number - 1, ind_2]) +
                                                     'GHz$', ha="center",
                                                     transform=plt.subplot(grid[ind_1, ind_2]).transAxes)

                if ind_2 != 0:
                    plt.subplot(grid[ind_1, ind_2]).axes.get_yaxis().set_visible(False)
                if ind_1 != 5:
                    plt.subplot(grid[ind_1, ind_2]).axes.get_xaxis().set_visible(False)
        plt.show()


    # show_modes_grid()
    #fmr(40, 0, 91, 1)
    make_final_plot()
    #cross_section(40)

