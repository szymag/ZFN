import matplotlib.pyplot as plt
import numpy as np


class Plot:
    def __init__(self, number_of_disp_branch, angle, name_of_output_file=None):
        self.name_of_file = name_of_output_file
        self. number_of_disp_branch = number_of_disp_branch + 1
        self.angle = angle

    def dispersion_relation(self, input_data):
        loaded_data = np.transpose(np.loadtxt(input_data))
        assert len(loaded_data.shape) == 2, 'loaded file must be two dimensional'
        for i in range(self.number_of_disp_branch):
            plt.plot(loaded_data[0] / 10e7, loaded_data[1 + i] / 10e8, color='r')
        #plt.ylim([2, 15])
        plt.xlabel(r'wektor odwrotny k $[nm^{-1}]$')
        plt.ylabel('częstotliwość [GHz]')
        self.show_or_save_plot()

    def fmr_freq_function_of_magnetic_field(self, begin_of_name_file,
                                            start_number, end_number, scaling_factor_x_axis=1):
        f = plt.figure()
        ax = f.add_subplot(111)
        plt.tight_layout(pad=6.5, w_pad=5.5, h_pad=5.0)
        x_axis = np.arange(start_number, end_number + 1) / scaling_factor_x_axis
        data = self.load_and_join_frequency_for_diff_field(begin_of_name_file, start_number, end_number)
        np.savetxt('fmr_'+str(self.angle)+'.txt', np.transpose(np.concatenate(([x_axis],data[ 1:]), axis=0 )))
        for i in range(1, self.number_of_disp_branch):
            if i%2 != 0:
                ax.plot(x_axis, data[i])
            else:
                ax.plot(x_axis, data[i], ls='--')
        ax.set_ylabel('frequency (GHz)', fontsize=22)
        ax.set_xlabel(r'external field $H_{0}$ (T)', fontsize=22)
        ax.set_title(r'FMR ($k=0$), $'+str(self.angle)+'^{\circ}$ ', fontsize=22)
        #ax.set_ylim([0, 20])
        #ax.set_xlim([0, 0.4])
        ax.xaxis.grid()
        ax.yaxis.grid()
        #self.create_legend_for_fmr(ax)
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
            plt.clf()
            plt.cla()
        else:
            plt.show()
            return 'wrong argument was puted'


if __name__ == "__main__":

    Plot(13, 0).dispersion_relation('dys_090.dat')
    #for i in np.array([0, 8, 16, 24, 32, 40]):
    #    Plot(10, 'dys_' + str(i)).dispersion_relation('fmr' + str(i) + '.dat')
    #for i in range(0, 92, 2):
        #Plot(3, i, 'fmr_'+str(i)+'deg').fmr_freq_function_of_magnetic_field(str(i) + '_', 1, 200, 1000)
    #tmp = np.loadtxt('freq_vs_angle_t900.dat')
    #print(len(np.arange(0, len(tmp[0,:]))))
    #print(len(tmp[0,:]))
    #for i in range(0, 50):
    #    plt.plot(np.arange(0, len(tmp[0,:])), tmp[i,:])
    #plt.ylim([1.3e9, 6e9])
    #plt.ylabel('frequency (GHz)', fontsize=22)
    #plt.xlabel('angle', fontsize=22)
    #plt.savefig('mode_freq_vs_angle_t900.eps')
