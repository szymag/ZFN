import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
from math import radians, sin
from src.eig_problem.ParametryMaterialowe import ParametryMaterialowe
from src.eig_problem.WektorySieciOdwrotnej import WektorySieciOdwrotnej


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
        magnetization = self.spatial_distribution_dynamic_magnetization(500, self.mode_number)
        elementary_cell = self.elementary_cell_reconstruction(500)
        fig = plt.figure()
        ax = fig.add_subplot(311)
        #ax2 = ax.twinx()
        ax.plot(magnetization[0], abs(magnetization[1])**2, '-', label=r'$\left|\mathbf{m}\right|^{2}$')
        #ax2.plot(magnetization[0], np.arctan2(magnetization[1].imag, magnetization[1].real), '-', label='phase', color="red")
        #ax2.plot(elementary_cell[0], elementary_cell[1], '-', label=r'$M_{s}$', color="green", linewidth=3)
        #ax.legend(loc=(0, .15), frameon=False)
        #ax2.legend(loc=(0, .05), frameon=False)
        #ax.grid()
        #ax.set_xlabel("elementary cell [nm]")
        #ax.set_ylabel(r"Intensity")
        #ax2.set_ylabel(r"Phase")
        #ax2.set_ylim(-np.pi, np.pi)
        #ax.set_ylim(0, 2)
        #ax.set_title(r'$angle = ' + str(self.angle) + '^{\circ}$, $mode = ' + str(self.mode_number + 1) + ', $'
        #             + '$H = 0.05T$', fontsize=22)

        #ax.set_ticks([])
        ax.set_xlabel('Position')
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
        x = np.linspace(0, self.lattice_const, grid)
        return sin(radians(2*self.angle))*np.trapz(tmp[1])

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
                Profile1D(mode, '/home/szymon/python/ZFN/src/eig_problem/0'+str(modulation)+'ni_' + str(int(i[1])) + '.dat',
                          'deg' + str(i[1]) + '_mode' + str(mode), angle=i[1]).generate_plot()
        else:
            for i in enumerate(np.arange(start, stop, step)):
                for j in range(1, 6):
                    Profile1D(j, '/home/szymon/python/ZFN/src/eig_problem/0'+str(modulation)+'ni_' + str(int(i[1])) + '.dat',
                              'deg' + str(i[1]) + '_mode' + str(j), angle=i[1]).generate_plot()


    def fmr(modulation, start, stop, step):
        fmr = np.zeros((50, int((stop - start) / step)))
        for i in enumerate(np.arange(start, stop, step)):
            fmr[:,i[0]] = Profile1D(1, '/home/szymon/python/ZFN/src/eig_problem/0'+str(modulation)+'ni_' + str(i[1]) + '.dat', None,
                  angle=i[1]).fmr_intensity_order()
        np.savetxt('fmr'+str(modulation)+'.dat', fmr)

    def fmr_intensity_map(modulation, start, stop, step):
        y = np.loadtxt('/home/szymon/python/ZFN/src/eig_problem/freq_vs_angle_0'+str(modulation)+'ni.dat') / 1e9
        if int((stop - start) / step) != y.shape[1]:
            raise ValueError('Current settings is different from ZagadnienieWlasne')
        x = np.arange(start, stop, step)
        try:
            z = np.loadtxt('fmr' + str(modulation) + '.dat')
            if z.shape[1] != y.shape[1]:
                fmr(modulation, start, stop, step)
                z = np.loadtxt('fmr' + str(modulation) + '.dat')
        except:
            fmr(modulation, start, stop, step)
            z = np.loadtxt('fmr' + str(modulation) + '.dat')
        max = np.max(z)
        min = np.min(z)
        for i in range(0, 13):
            plt.scatter(x, y[i], c=z[i,:], s=30, edgecolors='', vmin=min, vmax=max, alpha=0.9, cmap=plt.cm.binary)
        plt.ylim([4.8, 7.4])
        #plt.locator_params(nbins=30)
        #plt.xticks([])
        #plt.yticks([])
        #plt.show()

        plt.xlabel(r'angle [$^{\circ}$]', fontsize=20)
        plt.ylabel('frequency [GHz]', fontsize=20)
        a = plt.colorbar(ticks=None)
        a.set_ticks([])
        a.set_label('Arbitrary Unit [a.u.]')
        plt.savefig('fmr'+str(modulation)+'.svg')
        plt.close()
        plt.cla()

    def cross_section(modulation):
        fmr_map = np.loadtxt('fmr' + str(modulation) + '.dat')
        mode_freq = np.loadtxt('/home/szymon/python/ZFN/src/eig_problem/freq_vs_angle_0'+str(modulation)+'ni.dat')
        assert fmr_map.shape == mode_freq.shape
        points = np.zeros(fmr_map.shape[1], dtype=complex)
        for i in enumerate(np.arange(fmr_map.shape[1])):
            modes = np.where(np.abs(mode_freq[:, i[1]] - 5.196e9) < 1e9)[0]
            for j in modes:
                points[i[0]] += Profile1D(1, '/home/szymon/python/ZFN/src/eig_problem/0' + str(modulation) + 'ni_' + str(i[1]) + '.dat', None,
                      angle=i[1]).acoustic_cross_section(500, j) * mode_weight(mode_freq[j, i[1]])
        plt.plot(np.arange(fmr_map.shape[1]), np.abs(points))
        plt.savefig('cross_section' + str(modulation) + '.svg')
        plt.close()
        plt.cla()

    def mode_weight(position):
        peak_position = 5.196e9
        fwhm = 0.2e9
        return 4e15 / ((peak_position - position)**2 + (fwhm/2)**2)

    #for i in range(3,9):
    #    cross_section(i)
    #fmr(8, 0, 91, 1)
    #fmr_intensity_map(8, 0, 91, 1)
    #modes(8, 1, 30, 1, True)
    modes(8, 70, 90, 5, True)
    #fmr(8, 0, 91, 1)
    #cross_section(8)
