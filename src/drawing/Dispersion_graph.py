import matplotlib.pyplot as plt
import numpy as np

def plot_dispersion(file):
    print(file)
    file1 = np.transpose(np.loadtxt(file))
    for i in range(18):
        plt.plot(file1[0], file1[1 + i], color='r')
    #plt.xlabel('wektor falowy q [m^-1]')
    #plt.ylabel('f [Hz]')
    plt.title(file[1])
    #plt.ylim([0.7e10, 1.2e10])
    #plt.savefig('fig_' + str(file[0]) + '.svg')
    plt.show()
   # plt.close()
   # plt.clf()

#plot_dispersion('dys.txt')
#for i in enumerate(file):
#    plot_dispersion(i)

def magnetization_profile(down_level):
    x_axis = np.linspace(0, 2*np.pi, 1000)
    y_axis = (np.cos(x_axis) + 1) / 2 * 0.43 +0.57
    plt.ylim([0,1])
    plt.plot(x_axis, y_axis)
    plt.xlabel('normalized position')
    plt.ylabel('magnetization a.u.')
    plt.xticks([], [])
    plt.savefig('0.57.svg')

#magnetization_profile(1)

def frequency_k_0_versus_magnetic_field():
    x_axis = []
    y_axis_1 = []
    y_axis_2 = []
    y_axis_3 = []
    y_axis_4 = []
    y_axis_5 = []
    y_axis_6 = []
    y_axis_7 = []
    y_axis_8 = []
    y_axis_9 = []

    f = plt.figure()
    ax = f.add_subplot(111)
    plt.tight_layout(pad=6.5, w_pad=5.5, h_pad=5.0)
    for i in range(1, 400):
        x_axis.append(i)
        y_axis_1.append((np.loadtxt('0.5cdys_' + str(i) + '.dat') / 10e8)[1])
        y_axis_2.append((np.loadtxt('0.5cdys_' + str(i) + '.dat') / 10e8)[2])
        y_axis_3.append((np.loadtxt('0.5cdys_' + str(i) + '.dat') / 10e8)[3])
        y_axis_4.append((np.loadtxt('0.5cdys_' + str(i) + '.dat') / 10e8)[4])
        y_axis_5.append((np.loadtxt('0.5cdys_' + str(i) + '.dat') / 10e8)[5])
        y_axis_6.append((np.loadtxt('0.5cdys_' + str(i) + '.dat') / 10e8)[6])
        y_axis_7.append((np.loadtxt('0.5cdys_' + str(i) + '.dat') / 10e8)[7])
        y_axis_8.append((np.loadtxt('0.5cdys_' + str(i) + '.dat') / 10e8)[8])
        y_axis_9.append((np.loadtxt('0.5cdys_' + str(i) + '.dat') / 10e8)[9])
    ax.plot(np.array(x_axis)/1000, y_axis_1, linewidth=1.5, color='black')
    ax.plot(np.array(x_axis)/1000, y_axis_2, linewidth=1.5, color='black')
    ax.plot(np.array(x_axis)/1000, y_axis_3, linewidth=1.5, color='black')
    ax.plot(np.array(x_axis)/1000, y_axis_4, linewidth=1.5, color='black')
    ax.plot(np.array(x_axis)/1000, y_axis_5, linewidth=1.5, color='black')
    ax.plot(np.array(x_axis)/1000, y_axis_6, linewidth=1.5, color='black')
    ax.plot(np.array(x_axis)/1000, y_axis_7, linewidth=1.5, color='black')
    ax.plot(np.array(x_axis)/1000, y_axis_8, linewidth=1.5, color='black')
    ax.plot(np.array(x_axis)/1000, y_axis_9, linewidth=1.5, color='black')

    ax.set_ylabel('frequency (GHz)', fontsize=26)
    #ax.set_xticklabels((r'$\mathrm{0}$', ' ', ' ', ' ', ' ', '0.5', ' ', ' ', ' ', ' ', '1'))
    ax.set_xlabel('field (T)', fontsize=26)
    #3ax.set_xlabel('external field 'r"$\-H_{0}$"' (T)', fontsize=26)
    #ax.set_title('resonance frequency (k=0) as a function of field in backward-volume configuration')
    #ax.set_xticks(np.arange(0, 1.1, 0.1), minor=False)
    ax.xaxis.grid()
    ax.yaxis.grid()
    plt.show()
    #plt.savefig('bac_vol.svg')

frequency_k_0_versus_magnetic_field()

