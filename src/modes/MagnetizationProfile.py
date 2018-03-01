import matplotlib.pyplot as plt
import numpy as np

from src.eig_problem.InputParameter import InputParameter
from src.eig_problem.ReciprocalVector import ReciprocalVector


class Profile2D:
    def __init__(self, mode_number, load_data):
        # TODO: Eigvalue problem return fourier coefficients. Name of variables wrongly suggest that they're vectors.
        self.fourier_coefficient = np.loadtxt(load_data).view(complex)
        self.lattice_const_x = InputParameter.a
        self.lattice_const_y = InputParameter.b
        self.mode_number = mode_number - 1

    def generate_plot(self):
        pass

    def spatial_distribution_dynamic_magnetization(self, grid):
        mode = self.fourier_coefficient[self.mode_number, :81]
        self.inverse_discrete_fourier_transform(mode, [0, 0])


    def inverse_discrete_fourier_transform(self, data, vector_position):
        reciprocal_vectors = 2 * np.pi * ReciprocalVector(max(data.shape)).lista_wektorow2d('min')
        abs(np.sum(data * np.exp(1j * np.prod(reciprocal_vectors * vector_position, axis=1))))


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

    def generate_plot(self):
        magnetization = self.spatial_distribution_dynamic_magnetization(500)
        elementary_cell = self.elementary_cell_reconstruction(500)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(magnetization[0], abs(magnetization[1]) ** 2, '-', label=r'$\left|\mathbf{m}\right|^{2}$')
        #ax.plot(manetization[0], np.arctan2(magnetization[1].imag, magnetization[1].real),
        #        '-', label='phase', color="red")
        ax2 = ax.twinx()
        ax2.plot(elementary_cell[0], elementary_cell[1], '-', label=r'$M_{s}$', color="green", linewidth=3)
        ax.legend(loc=(0, .1), frameon=False)
        ax2.legend(loc=(0, .05), frameon=False)
        ax.grid()
        ax.set_xlabel("elementary cell [nm]")
        ax.set_ylabel(r"Intensity")
        ax2.set_ylabel(r"Magnetization saturation $M_{s,Ni}\left(x\right)/M_{s,Ni}$")
        ax2.set_ylim(0, 1)
        ax.set_ylim(0, 2)
        ax.set_title(r'$angle = ' + str(self.angle) + '^{\circ}$, $mode = ' + str(self.mode_number + 1) + ', $'
                     + '$H = 0.05T$', fontsize=22)
        self.output_plot()

    def output_plot(self):
        if self.name_of_file is None:
            plt.show()
        elif type(self.name_of_file) == str:
            plt.savefig(self.name_of_file + '_' + '.png')
            plt.clf()
            plt.close()
        else:
            plt.show()
            return 'wrong argument was puted'

    def save_to_file(self):

        to_file = self.spatial_distribution_dynamic_magnetization(500)[1].real, \
                  self.spatial_distribution_dynamic_magnetization(500)[1].imag

        np.savetxt(self.name_of_file + '.txt', np.transpose(to_file))

    def spatial_distribution_dynamic_magnetization(self, grid):
        mode = self.eig_vectors[self.mode_number, 0:50]
        x = np.linspace(-self.lattice_const, 0, grid)
        tmp = np.zeros(grid, dtype=complex)
        for ind in enumerate(x):
            tmp[ind[0]] = self.inverse_discrete_fourier_transform(mode, ind[1])
        return x * 10 ** 9, tmp

    def elementary_cell_reconstruction(self, grid):
        coefficient = np.transpose(np.loadtxt('c_coef_100.txt').view(complex))
        x = np.linspace(-self.lattice_const, 0, grid)
        tmp = np.zeros(grid)
        for ind in enumerate(x):
            tmp[ind[0]] = abs(self.inverse_discrete_fourier_transform(coefficient, ind[1]))
        return x * 10 ** 9, tmp / (10 / 7) + 0.3

    def inverse_discrete_fourier_transform(self, data, vector_position):
        reciprocal_vectors = np.array(2 * np.pi * ReciprocalVector(max(data.shape)).lista_wektorow1d('min')
                                      / self.lattice_const)
        return np.sum(data * np.exp(1j * reciprocal_vectors * vector_position))


if __name__ == "__main__":
    #for i in range(0, 93, 3):
    #    for j in range(1, 2):
    #        Profile1D(j, 'vec_' + str(i) + '.dat', 'mode' + str(j) + '_deg' + str(i), angle=i).generate_plot()

    Profile2D(1, 'tst.vec').spatial_distribution_dynamic_magnetization(50)
