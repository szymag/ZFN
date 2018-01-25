import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
from math import radians, sin
from src.eig_problem.ParametryMaterialowe import ParametryMaterialowe
from src.eig_problem.WektorySieciOdwrotnej import WektorySieciOdwrotnej
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

class Profile1D:
    def __init__(self, mode_number, load_data, name_of_file, **kwargs):
        self.eig_vectors = np.loadtxt(load_data).view(complex)
        self.lattice_const = ParametryMaterialowe.a
        self.mode_number = mode_number - 1
        self.name_of_file = name_of_file

        if 'angle' in kwargs:
            self.angle = kwargs['angle']
        if 'field' in kwargs:
            self.field = kwargs['field']

    def generate_plot(self, ax, color_index):
        colors = ['C0', 'C1']
        x, magnetization = self.spatial_distribution_dynamic_magnetization(500, self.mode_number)
        phase = np.abs(np.arctan2(magnetization.imag, magnetization.real))
        parameter = np.array([0 if i < 0.5*np.pi else 1 for i in phase])
        magnetization1 = np.ma.masked_where(parameter != 0, abs(magnetization)**2)
        magnetization2 = np.ma.masked_where(parameter == 0, abs(magnetization)**2)
        ax.plot(x, magnetization1, colors[color_index] + '-',
                x, magnetization2, colors[color_index] + '--')

        ax.set_ylim(-0.05, 3)
        plt.setp(ax, xticks=[-1100, 0, 1100], xticklabels=[r'$-\Lambda$', 0, r'$\Lambda$'],
                 yticks=[0])
        ax.set_xlabel('Position')
        ax.set_ylabel(r'$\left|m_{out}\right|$')
        #self.output_plot()

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
        mode = self.eig_vectors[mode_number, 0:self.eig_vectors.shape[1]//2]
        x = np.linspace(-self.lattice_const, self.lattice_const, grid)
        tmp = np.zeros(grid, dtype=complex)
        for ind in enumerate(x):
            tmp[ind[0]] = self.inverse_discrete_fourier_transform(mode, ind[1])
        return x* 10 ** 9, tmp

    def elementary_cell_reconstruction(self, grid):
        coefficient = np.transpose(np.loadtxt('c_coef_100.txt').view(complex))
        x = np.linspace(-self.lattice_const, 0, grid)
        tmp = np.zeros(grid)
        for ind in enumerate(x):
            tmp[ind[0]] = abs(self.inverse_discrete_fourier_transform(coefficient, ind[1]))
        return x* 10 ** 9, tmp / (10 / 7) + 0.3

    def inverse_discrete_fourier_transform(self, data, vector_position):
        reciprocal_vectors = np.array(2 * np.pi * WektorySieciOdwrotnej(max(data.shape)).lista_wektorow1d('min') / self.lattice_const)
        return np.sum(data * np.exp(1j * reciprocal_vectors * vector_position))

    def fmr_intensity(self, mode):
        return np.abs(np.trapz(mode))**2 / np.trapz(np.abs(mode)**2)

    def fmr_intensity_order(self):
        fmr = [self.fmr_intensity(self.spatial_distribution_dynamic_magnetization(500, i)[1]) for i in range(50)]
        #np.savetxt(self.name_of_file, np.array(fmr))
        return np.array(fmr)

    def acoustic_cross_section(self, grid, data=None):
        tmp = self.spatial_distribution_dynamic_magnetization(grid, data)
        x = np.linspace(-self.lattice_const, self.lattice_const, grid)
        return sin(radians(2*self.angle))*np.trapz(tmp[1] * np.cos(x*2*np.pi / self.lattice_const))

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
                    Profile1D(mode[1], '/home/szymag/python/ZFN/src/eig_problem/' + str(int(i)) + '.dat',
                          'deg' + str(i) + '_mode' + str(mode[1]), angle=i).generate_plot(ax, mode[0])
            #x = np.linspace(-1100, 1100, 500)
            #ax.plot(x, np.cos(2*np.pi*x / 1100)+1, ls='--')
        else:
            for i in enumerate(np.arange(start, stop, step)):
                for j in range(5, 8):
                    Profile1D(j, '/home/szymag/python/ZFN/src/eig_problem/'+ str(int(i[1])) + '.dat',
                              'deg' + str(i[1]) + '_mode' + str(j), angle=i[1]).generate_plot(ax, i[0])


    def fmr(modulation, start, stop, step):
        fmr = np.zeros((50, int((stop - start) / step)))
        for i in enumerate(np.arange(start, stop, step)):
            fmr[:, i[0]] = Profile1D(1, '/home/szymag/python/ZFN/src/eig_problem/'+ str(i[1]) + '.dat', None,
                  angle=i[1]).fmr_intensity_order()
        np.savetxt('fmr'+str(modulation)+'.dat', fmr)


    def fmr_1(modulation, field, start, stop, step):
        return Profile1D(1, '/home/szymag/python/ZFN/src/eig_problem/' + str(field)
                                     + '_' + str(modulation) + '.dat', None).fmr_intensity_order()


    def fmr_intensity_map(modulation, start, stop, step, ax):
        #mateusz_plot(ax)
        if type(step) is float:
            y = np.loadtxt('/home/szymag/python/ZFN/src/eig_problem/densefreq_vs_angle_'+str(modulation)+'.dat') / 1e9
        else:
            y = np.loadtxt('/home/szymag/python/ZFN/src/eig_problem/freq_vs_angle_'+str(modulation)+'.dat') / 1e9

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
            tmp = ax.scatter(np.zeros(drawing_mode_count)+i[1],
                        y[0:drawing_mode_count, i[0]][order],
                        c=z[0:drawing_mode_count, i[0]][order],
                        s=15, edgecolors='', vmin=min, vmax=max,
                        cmap=plt.cm.Wistia, alpha=0.7, norm=MidpointNormalize(midpoint=30))
        ax.set_ylim([4.2, 10])
        #plt.ylim([5.06, 5.3])
        #ax.set_ylim.locator_params(nbins=11)
        #plt.xticks([13, 19.9], [13, 20])
        #plt.yticks([5.06, 5.3], [5.06, 5.3])
        ax.set_xlabel('Angle (Deg)')
        ax.set_ylabel('Frequency (GHz)')
        plt.colorbar(tmp, orientation="horizontal", ax=ax, ticks=[])
        #a.set_label('FMR Intensity (a.u.)')
        #if type(step) is float:
        #    plt.savefig('densefmr' + str(modulation) + '.svg')
        #elif type(step) is int:
        #    plt.savefig('fmr'+str(modulation)+'.svg')
        #plt.close()
        #plt.cla()


    def cross_section(modulation):
        fmr_map = np.loadtxt('fmr' + str(modulation) + '.dat')
        mode_freq = np.loadtxt('/home/szymag/python/ZFN/src/eig_problem/freq_vs_angle_'+str(modulation)+'.dat')
        assert fmr_map.shape == mode_freq.shape
        points = np.zeros(fmr_map.shape[1], dtype=complex)
        for i in enumerate(np.arange(fmr_map.shape[1])):
            mode = np.argsort(fmr_map[:, i[1]])[-1]
            points[i[0]] += Profile1D(mode, '/home/szymag/python/ZFN/src/eig_problem/' + str(i[1]) + '.dat',
                                              None, angle=i[1]).acoustic_cross_section(500, mode) * mode_weight(mode_freq[mode, i[1]])
        fig = plt.figure()
        ax = fig.add_subplot(211)
        ax.plot(np.arange(0, 91), np.abs(points))
        ax.set_xlabel('Angle (Deg)')
        ax.set_ylabel('Arbitrary Unit (a.u.)')
        ax.set_yticks([], [])
        plt.savefig('cross.section.png')


    def mode_weight(position):
        peak_position = 4.8e9
        fwhm = 0.3e9
        return 8e15 / ((peak_position - position)**2 + (fwhm/2)**2)

    def mateusz_plot(ax):
        mateusz = np.loadtxt('szymon_fmr.txt', delimiter=',')[:, 300:600]
        freq_mat = np.loadtxt('szymon_freq.txt')[300:600] / 1e9

        x = np.linspace(0, 90, 45)
        X, Y = np.meshgrid(freq_mat, x)
        ax.pcolor(Y, X, mateusz, cmap=plt.cm.binary)
        #plt.locator_params(nbins=11)
        plt.ylim([4.4, 7])
        #plt.xticks([], [])
        #plt.yticks([], [])
        #a = plt.colorbar(ticks=[])

        #a.set_label('FMR Intensity (a.u.)')
        #plt.savefig('mateusz.svg')
        #plt.show()

    #cross_section(155)

    #fmr(40, 0, 91, 1)
    #fmr_intensity_map(155, 0, 91, 1)
    #modes('heat', 0, 90, 1, True)
    #modes(155, 0, 70, 2, True)
    #modes(8, 16, 17, 1, False)
    #fmr(8, 0, 91, 1)
    #modulations = [948, 935, 922, 915]
    #for i in modulations:
    #    plt.plot(*cross_section(i), label='Max Ms = 0.' + str(i) + 'Ni')
    #plt.legend()
    #plt.savefig('cross_section' '.png')
    #mateusz_plot()
    import matplotlib.gridspec as gridspec
    from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
    from mpl_toolkits.axes_grid1.inset_locator import mark_inset

    def get_y_coordinates(angles, modulation):
        fmr_map = np.loadtxt('fmr' + str(modulation) + '.dat')
        mode_freq = np.loadtxt('/home/szymag/python/ZFN/src/eig_problem/freq_vs_angle_' + str(modulation) + '.dat')
        return [mode_freq[np.argsort(fmr_map[:, i])[-2:], i] / 1e9 for i in angles]

    def make_final_plot():
        angles = [0, 10, 20, 30, 40, 60]
        y_coordinates = get_y_coordinates(angles, 40)
        labels_1 = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)']
        labels_2 = ['(g)', '(h)', '(i)', '(j)', '(k)', '(l)']

        outer = gridspec.GridSpec(2, 1, height_ratios=[4, 6])

        gs1 = gridspec.GridSpecFromSubplotSpec(3, 2, subplot_spec = outer[0],
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

        gs2 = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=outer[1])
        ax7 = plt.subplot(gs2[0, :])
        fmr_intensity_map(40, 0, 91, 1, ax7)
        for i, j, k in zip(angles, labels_1, y_coordinates):
            ax7.text(i, k[1], j, ha="center")
        for i, j, k in zip(angles, labels_2, y_coordinates):
            ax7.text(i, k[0], j, ha="center")
        plt.show()


    def show_modes_grid():
        fields = [0.03, 0.035, 0.04, 0.045, 0.05, 0.055, 0.06, 0.065]
        modulations = [793, 837, 877, 894, 903, 909]

        grid = gridspec.GridSpec(6, 8, hspace=0.0, wspace=0.0)
        f = plt.figure()
        f.suptitle('Angle=60 Deg')
        #ax1.title.set_text('sdf')
        for ind_1, modulation in enumerate(modulations):
            frequencies = np.around(np.loadtxt('/home/szymag/python/ZFN/src/eig_problem/freq_vs_angle_' + str(modulation) + '.dat') / 1e9, decimals=2)
            for ind_2, field in enumerate(fields):
                mode_number = np.argsort(fmr_1(modulation, field, 0, 1, 1))[-1] + 1
                Profile1D(mode_number, '/home/szymag/python/ZFN/src/eig_problem/'
                          + str(field) + '_' + str(modulation) + '.dat', 'dummy',
                          field=field, angle=60).generate_plot(plt.subplot(grid[ind_1, ind_2]), 0)

                plt.subplot(grid[ind_1, ind_2]).text(0.2, 0.1, r'$H_{0}=' + str(field) +'$', ha="center",
                        transform=plt.subplot(grid[ind_1, ind_2]).transAxes)
                plt.subplot(grid[ind_1, ind_2]).text(0.2, 0.8, r'$M_{min}=0.' + str(modulation)+'$', ha="center",
                        transform=plt.subplot(grid[ind_1, ind_2]).transAxes)
                plt.subplot(grid[ind_1, ind_2]).text(0.5, 0.5, r'$freq=' +
                                                     str(frequencies[mode_number-1, ind_2]) +
                                                     'GHz$', ha="center",
                        transform=plt.subplot(grid[ind_1, ind_2]).transAxes)

                if ind_2 != 0:
                    plt.subplot(grid[ind_1, ind_2]).axes.get_yaxis().set_visible(False)
                if ind_1 != 5:
                    plt.subplot(grid[ind_1, ind_2]).axes.get_xaxis().set_visible(False)
        plt.show()

    #show_modes_grid()
    make_final_plot()