import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
from math import radians, sin
from src.eig_problem.ParametryMaterialowe import ParametryMaterialowe
from src.eig_problem.WektorySieciOdwrotnej import WektorySieciOdwrotnej
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

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

    def generate_plot(self):
        magnetization1 = self.spatial_distribution_dynamic_magnetization(500, self.mode_number)

        elementary_cell = self.elementary_cell_reconstruction(500)
        fig = plt.figure()
        ax = fig.add_subplot(311)
        #ax2 = ax.twinx()
        ax.plot(magnetization1[0], abs(magnetization1[1])**2, '-', label=r'$\left|\mathbf{m}\right|^{2}$', color='black')

        #ax2.plot(magnetization[0], np.arctan2(magnetization[1].imag, magnetization[1].real), '-', label='phase', color="red")
        #ax2.plot(elementary_cell[0], elementary_cell[1], '-', label=r'$M_{s}$', color="green", linewidth=3)
        #ax.legend(loc=(0, .15), frameon=False)
        #ax2.legend(loc=(0, .05), frameon=False)
        #ax.grid()
        #ax.set_xlabel("elementary cell [nm]")
        #ax.set_ylabel(r"Intensity")
        #ax2.set_ylabel(r"Phase")
        #ax2.set_ylim(-np.pi, np.pi)
        ax.set_ylim(0, 3)
        #ax.set_title(r'$angle = ' + str(self.angle) + '^{\circ}$, $mode = ' + str(self.mode_number + 1) + ', $'
        #             + '$H = 0.05T$', fontsize=22)
        #plt.locator_params(nbins=3)
        plt.xticks([-1100, 0, 1100], [r'$-\Lambda$', 0, r'$\Lambda$'])
        #ax.set_ticks([])
        #plt.yticks([0], [0])
        ax.set_xlabel('Position')
        ax.set_ylabel('$|m|$ (a.u.)')
        self.output_plot()

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
        x = np.linspace(0, 2*self.lattice_const, grid)
        return sin(radians(2*self.angle))*np.trapz(tmp[1] * np.cos(x*2*np.pi / self.lattice_const))

if __name__ == "__main__":

    def modes(modulation, start, stop, step, only_max_fmr=True):
        if only_max_fmr:
            try:
                fmr_map = np.loadtxt('fmr' + str(modulation) + '.dat')
            except:
                fmr(modulation, start, stop, step)
                fmr_map = np.loadtxt('fmr' + str(modulation) + '.dat')

            for i in enumerate(np.arange(start, stop, step)):
                mode = np.argmax(fmr_map[:, i[1]]) + 1
                Profile1D(mode, '/home/szymag/python/ZFN/src/eig_problem/'+ str(int(i[1])) + '.dat',
                          'deg' + str(i[1]) + '_mode' + str(mode), angle=i[1]).generate_plot()
        else:
            for i in enumerate(np.arange(start, stop, step)):
                for j in range(5, 6):
                    Profile1D(j, '/home/szymag/python/ZFN/src/eig_problem/'+ str(int(i[1])) + '.dat',
                              'deg' + str(i[1]) + '_mode' + str(j), angle=i[1]).generate_plot()


    def fmr(modulation, start, stop, step):
        fmr = np.zeros((50, int((stop - start) / step)))
        for i in enumerate(np.arange(start, stop, step)):
            fmr[:,i[0]] = Profile1D(1, '/home/szymag/python/ZFN/src/eig_problem/'+ str(i[1]) + '.dat', None,
                  angle=i[1]).fmr_intensity_order()
        np.savetxt('fmr'+str(modulation)+'.dat', fmr)

    def fmr_intensity_map(modulation, start, stop, step):
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
        drawing_mode_count = 20
        for i in enumerate(np.arange(start, stop, step)):
            order = z[0:drawing_mode_count,i[0]].argsort()
            plt.scatter(np.zeros(drawing_mode_count)+i[1],
                        y[0:drawing_mode_count,i[0]][order],
                        c=z[0:drawing_mode_count,i[0]][order],
                        s=15, edgecolors='', vmin=min, vmax=max,
                        cmap=plt.cm.Blues, alpha=0.7,norm=MidpointNormalize(midpoint=50))
        plt.ylim([4.4, 7])
        #plt.ylim([5.06, 5.3])
        plt.locator_params(nbins=11)
        #plt.xticks([13, 19.9], [13, 20])
        #plt.yticks([5.06, 5.3], [5.06, 5.3])
        plt.xlabel('Angle (Deg)', fontsize=20)
        plt.ylabel('Frequency (GHz)', fontsize=20)
        a = plt.colorbar(ticks=[])

        a.set_label('FMR Intensity (a.u.)')
        if type(step) is float:
            plt.savefig('densefmr' + str(modulation) + '.svg')
        elif type(step) is int:
            plt.savefig('fmr'+str(modulation)+'.svg')
        plt.close()
        plt.cla()


    def cross_section(modulation):
        fmr_map = np.loadtxt('fmr' + str(modulation) + '.dat')
        mode_freq = np.loadtxt('/home/szymag/python/ZFN/src/eig_problem/freq_vs_angle_'+str(modulation)+'.dat')
        assert fmr_map.shape == mode_freq.shape
        points = np.zeros(fmr_map.shape[1], dtype=complex)
        for i in enumerate(np.arange(fmr_map.shape[1])):
            mode = np.argsort(fmr_map[:, i[1]])[-2:]
            for j in mode:
                points[i[0]] = Profile1D(j, '/home/szymag/python/ZFN/src/eig_problem/' + str(i[1]) + '.dat',
                                              None, angle=i[1]).acoustic_cross_section(500, j)
        plt.plot(np.arange(0, 91), np.abs(points))
        plt.show()


    def mode_weight(position):
        peak_position = 4.8e9
        fwhm = 0.2e9
        return 4e15 / ((peak_position - position)**2 + (fwhm/2)**2)

    def mateusz_plot():
        mateusz = np.loadtxt('szymon_fmr.txt', delimiter=',')[:,:100] / 1e9
        freq_mat = np.loadtxt('szymon_freq.txt')[100] / 1e9

        x = np.linspace(0, 90, 45)
        y = np.linspace(0, freq_mat, 100)
        X, Y = np.meshgrid(y, x)
        plt.pcolor(Y, X, mateusz, cmap=plt.cm.binary)
        plt.locator_params(nbins=11)
        plt.ylim([4.4, 7])
        #plt.xticks([], [])
        #plt.yticks([], [])
        a = plt.colorbar(ticks=[])

        a.set_label('FMR Intensity (a.u.)')
        plt.savefig('mateusz.svg')

    #cross_section(8)
    #for i in range(3,9):
    #    cross_section(i)
    #fmr(8, 0, 91, 1)
    #fmr_intensity_map('heat', 0, 91, 1)
    #fmr_intensity_map(8, 0, 91, 1)
    #modes('heat', 0, 90, 1, True)
    modes(8, 0, 40, 5, True)
    #modes(8, 16, 17, 1, False)
    #fmr(8, 0, 91, 1)
    #modulations = [948, 935, 922, 915]
    #for i in modulations:
    #    plt.plot(*cross_section(i), label='Max Ms = 0.' + str(i) + 'Ni')
    #plt.legend()
    #plt.savefig('cross_section' '.png')
    mateusz_plot()