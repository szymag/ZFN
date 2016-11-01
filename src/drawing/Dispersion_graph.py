import matplotlib.pyplot as plt
import numpy as np

name1 = 'Ni0.5Ni_0.04.txt'
name2 = 'PyCo_c.txt'
name3 = 'NiNi_p.txt'
name4 = 'NiNi_c.txt'
name5 = 'Ni0.8Ni_p.txt'
name6 = 'Ni0.8Ni_c.txt'
name7 = 'Ni0.5Ni_p.txt'
name8 = 'Ni0.5Ni_c.txt'
#file = [name1, name2, name3, name4, name5, name6, name7, name8]
file = [name1]
def plot_dispersion(file):
    print(file)
    file1 = np.transpose(np.loadtxt(file[1]))
    for i in range(5):
        plt.plot(file1[0], file1[1 + i], color='r')
    plt.xlabel('wektor falowy q [m^-1]')
    plt.ylabel('f [Hz]')
    plt.title(file[1])
    #plt.ylim([0.7e10, 1.2e10])
    #plt.savefig('fig_' + str(file[0]) + '.svg')
    plt.show()
    plt.close()
    plt.clf()


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
    for i in np.arange(0.01, 0.20, 0.01):
        x_axis.append(i)
        y_axis_1.append(np.loadtxt('Ni0.5Ni_' + str(i) + '.txt')[1])
        y_axis_2.append(np.loadtxt('Ni0.5Ni_' + str(i) + '.txt')[2])
        y_axis_3.append(np.loadtxt('Ni0.5Ni_' + str(i) + '.txt')[3])
        y_axis_4.append(np.loadtxt('Ni0.5Ni_' + str(i) + '.txt')[4])
        y_axis_5.append(np.loadtxt('Ni0.5Ni_' + str(i) + '.txt')[5])
    plt.plot(x_axis, y_axis_1)
    plt.plot(x_axis, y_axis_2)
    plt.plot(x_axis, y_axis_3)
    plt.plot(x_axis, y_axis_4)
    plt.plot(x_axis, y_axis_5)
    plt.show()

frequency_k_0_versus_magnetic_field()

