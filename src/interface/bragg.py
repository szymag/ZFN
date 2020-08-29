import numpy as np
from scipy.signal import find_peaks
from src.fft_from_image.FFT import FFT
import matplotlib.pyplot as plt
import matplotlib
lattice_const = 34307e-9

# exchange constant
A = (1.1e-11*2)/2

# magnetization saturation
ms = (144*1.445e6 + 233*0.86e6)/377
# thickness
d = 30e-9

mu0 = 1.2566370614359173e-06
gamma = 176e9

lamb2 = (144*3.46e-17+233*2.80e-17)/377
# external magnetic field
H = 0.1


# wave vector
k = np.linspace(1e1,0.4e8, 400)

def dispersion(k):
    Nd=0.0
    omega0 = gamma*H
    omegaM = gamma*mu0*ms
    el_1 = omega0 + omegaM*lamb2*k*k + Nd*omegaM
    el_2 = omega0 + omegaM*lamb2*k*k + (1-Nd)*omegaM
    el_3 = omegaM**2/4*(1 - np.exp(-2*d*k))
    return np.sqrt(el_1*el_2+el_3) / 2 / np.pi / 1e9


# dys = np.loadtxt('./src/interface/Fib_14.dys') /1e9
dys = np.loadtxt('idos')/1e9


# phason_dys = np.loadtxt('/media/szymag/Dane/ZFN/JEMS2019/fib/F_5_' + str(0) + '.dys') /1e9
# coef = abs(FFT().wywolaj_fft1d('P', 10, 377))
coef = np.transpose(np.loadtxt('/media/szymag/Dane/ZFN/JEMS2019/periodic/P_5_0.fft'))

coef = coef[0] + coef[1]*1j
struct = abs(np.fft.ifft(coef))
coef = coef[len(coef)//2:]/max(coef)*10
coef = abs(coef[:])


x_coef = 2*np.pi*np.arange(0, len(coef))/lattice_const
x = np.linspace(0, 1e8, len(dys))


def get_peaks(height):
    peaks = find_peaks(coef, height)[0]
    for num, i in enumerate(peaks):
        print(str(num) + ' peak at ' + str(x_coef[i]) + ' with height ' + str(coef[i]))
    return peaks


def plot_peaks_idos():
    peaks = get_peaks(0.7)
    print(peaks)

    # plt.vlines(x_coef[peaks], 0, 30)
    # plt.scatter(x_coef[peaks], coef[peaks], color='r')
    plt.scatter(x, dys)
    # plt.plot(x_coef, coef)
    plt.show()


def plot_gaps():
    relative_height = np.zeros(100)
    relative_width = np.zeros(100)
    peaks = get_peaks(0.7)
    peak_num = input('Which one? ')
    peak = peaks[int(peak_num)]

    for i in range(100):
        phason_fft = abs(np.loadtxt('/media/szymag/Dane/ZFN/JEMS2019/fib/F_5_' + str(i) + '.fft').view(complex))
        phason_fft = phason_fft[len(phason_fft)//2:]/max(phason_fft)*10
        phason_fft = phason_fft[:len(dys)]
        phason_dys = np.loadtxt('/media/szymag/Dane/ZFN/JEMS2019/fib/F_5_' + str(i) + '.dys') /1e9
        # print(peaks)
        relative_height[i] = phason_fft[peak] / coef[peak]
        relative_width[i] = phason_dys[peak] - phason_dys[peak-1]
    plt.xlabel('Realtive value of Fourier components (%)')
    plt.ylabel('Gap width (GHz)')
    plt.scatter(relative_height, relative_width)
    plt.show()

def averaged_peaks(peak_pos):
    full_peak = abs(np.loadtxt('./f_coef_10*14.fft').view(complex))
    full_peak = full_peak[len(full_peak)//2:]/max(full_peak)
    full_peak = full_peak[:len(dys)]
    full_peak = full_peak[peak_pos]
    h = np.zeros(10)
    s = np.zeros(10)
    for j, phas in enumerate([5, 15, 25, 35, 45, 55, 65, 80, 100, 120]):
        height = np.zeros(100)
        for i in range(100):
            phason_fft = abs(np.loadtxt('/media/szymag/Dane/ZFN/JEMS2019/fib/F_'+str(phas)+'_' + str(i) + '.fft').view(complex))

            phason_fft = phason_fft[len(phason_fft)//2:]/max(phason_fft)
            phason_fft = phason_fft[:len(dys)]
            # plt.plot(phason_fft)
            # plt.show()
            height[i] = phason_fft[peak_pos]
        h[j] = np.mean(1-height/full_peak)
        s[j] = np.std(1-height/full_peak)
    return 100*h, 100*s

def averaged_gaps(gap_width, rang):

    g = np.zeros(10)
    s = np.zeros(10)
    for j, phas in enumerate([5, 15, 25, 35, 45, 55, 65, 80, 100, 120]):
        gap = np.zeros(100)
        for i in range(100):
            phason_dys = np.loadtxt('/media/szymag/Dane/ZFN/JEMS2019/fib/F_'+str(phas)+'_' + str(i) + '.dys')[rang[0]:rang[1]] /1e9

            # print(np.max(np.diff(phason_dys)))
            gap[i] = np.max(np.diff(phason_dys))
            # plt.scatter(np.arange(len(phason_dys)), phason_dys)
            # plt.show()

        g[j] = np.mean(1-gap/gap_width)
        s[j] = np.std(1-gap/gap_width)
    return 100*g, 100*s

def gap_evolution():
    plt.style.use('seaborn-colorblind')
    plt.rc('text', usetex=False)
    matplotlib.rcParams['font.family'] = "Liberation Sans"
    fig = plt.figure(figsize=(3.3,2.2))
    full_gap_1 = 0.52
    full_gap_2 = 3.22
    full_gap_3 = 1.62
    g, sy = averaged_gaps(full_gap_1, [80, 100])
    h, sx = averaged_peaks(89)
    print(g)
    plt.errorbar(h, g, xerr=sx, yerr=sy, fmt='o', label='Third widest gap')
    g, sy = averaged_gaps(full_gap_2, [120, 200])
    h, sx = averaged_peaks(144)
    print(g)
    plt.errorbar(h, g, xerr=sx, yerr=sy, fmt='o', label='First widest gap')
    g, sy = averaged_gaps(full_gap_3, [250, 350])
    h, sx = averaged_peaks(233)
    print(g)
    plt.errorbar(h, g, xerr=sx, yerr=sy, fmt='o', label='Second widest gap')

    plt.legend()
    plt.xlabel('Peak height reduction (%)')
    plt.ylabel('Gap width reduction (%)')

    # plt.show()
    plt.tight_layout()
    plt.savefig('width.svg')


def bragg_sub(series_count):
    plt.style.use('seaborn-colorblind')
    plt.rc('text', usetex=False)
    fig = plt.figure(figsize=(14,10))
    gs = fig.add_gridspec(2,2)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 0])
    # ax4 = fig.add_subplot(gs[2, 1])
    ax1.plot(k, dispersion(k))
    ax1.set_xlim([-0.3e7, 4e7])
    phason_dys = np.loadtxt('/media/szymag/Dane/ZFN/JEMS2019/fib/F_15_' + str(series_count) + '.dys') /1e9
    ax2.scatter(np.arange(len(phason_dys)), phason_dys, s=1)
    ax2.plot(dispersion(k))
    phason_fft = abs(np.loadtxt('/media/szymag/Dane/ZFN/JEMS2019/fib/F_15_' + str(series_count) + '.fft').view(complex))
    phason_fft = phason_fft[len(phason_fft)//2:]/max(phason_fft)
    phason_fft = phason_fft[:len(dys)]
    full_peak = abs(np.loadtxt('./f_coef_10*14.fft').view(complex))
    full_peak = full_peak[len(full_peak)//2:]/max(full_peak)
    full_peak = full_peak[:len(dys)]

    ax3.plot(2*np.pi*np.arange(0, len(phason_fft))/lattice_const, full_peak, alpha=0.5, ls='--', color='C0')
    ax3.plot(2*np.pi*np.arange(0, len(phason_fft))/lattice_const, phason_fft)
    ax3.set_xlim([-0.3e7,8e7])
    plt.tight_layout()
    # plt.show()


def count_full_gap(idos):
    basic_sep = (idos[1] - idos[0])/1.
    gaps = np.diff(idos)
    gaps = gaps[basic_sep < gaps]
    gaps = np.sum(gaps)
    return gaps/1e9


def plot_full_gap_evolution():
    plt.style.use('seaborn-colorblind')
    plt.rc('text', usetex=False)
    matplotlib.rcParams['font.family'] = "Liberation Sans"
    fig = plt.figure(figsize=(3.3,2.2))
    idos = np.loadtxt('./Fib_14.dys')
    reference_gap = count_full_gap(idos)
    print(reference_gap)
    gaps = np.zeros(10)
    for j, phas in enumerate([5, 15, 25, 35, 45, 55, 65, 80, 100, 120]):
        tmp = np.zeros(100)
        for i in range(100):
            idos = np.loadtxt('/media/szymag/Dane/ZFN/JEMS2019/fib/F_'+str(phas)+'_' + str(i) + '.dys')
            tmp[i] = count_full_gap(idos)
        gaps[j] = 1 - np.mean(tmp)/reference_gap
    plt.scatter([5, 15, 25, 35, 45, 55, 65, 80, 100, 120], gaps)
    plt.xlabel('Phasons')
    plt.ylabel('Gap width reduction (%)')
    plt.tight_layout()
    plt.savefig('gap_evolution.svg')


def hausdorff(idos, r):
    import matplotlib.patches as patches
    idos = idos/1e9
    # r *= 1e9

    count_in = np.zeros(len(idos))
    for idx in range(1, len(idos)-1):
        dist_right = idos[idx+1] - idos[idx]
        dist_left = idos[idx] - idos[idx-1]
        if dist_right/dist_left > 40:
            # print(dist_right/dist_left)
            count_in[idx] = 1
        if dist_left/dist_right > 20:
            count_in[idx] = 1
    count_in[0] = 1
    count_in[-1] = 1
    pass_bang_edges = np.argwhere(count_in)

    fig,ax = plt.subplots(1)
    count_hit = 0
    for idx in range(0, len(pass_bang_edges), 2):
        pass_left = idos[pass_bang_edges[idx]][0]
        pass_right = idos[pass_bang_edges[idx+1]][0]
        w = pass_right - pass_left
        h = 400

        n1 = (pass_left // r)*r + r
        n2 = (pass_right // r)*r + r

        if pass_left < n1 and pass_right <= n1:
            # print(pass_left, pass_right, 1)
            count_hit += 0
        elif pass_left < n1 and pass_right < n1+r:
            count_hit += 2
            print('asas')
            # print(pass_left, pass_right, 2)
        elif pass_left < n1 and pass_right > (n1+r):
            count_hit += int(1 + (n2-n1)//r)
            # print(pass_left, pass_right, int(1 + (n2-n1)//r))
        else:
            print('Someting wrong' + str(pass_left))
        plt.axvline(pass_left, 0, 400)
        plt.axvline(pass_right, 0, 400)
        ax.add_patch(patches.Rectangle((pass_left, 0), w, h, alpha=0.3))

    plt.step(idos, np.arange(len(idos)))
    plt.show()
    return count_hit


# idos = np.loadtxt('./Fib_14.dys')
# idos = np.loadtxt('./Periodic.dys')
# idos = np.loadtxt('/media/szymag/Dane/ZFN/JEMS2019/fib/F_'+str(120)+'_' + str(0) + '.dys')
# idos = np.loadtxt('/media/szymag/Dane/ZFN/JEMS2019/periodic/P_'+str(5)+'_' + str(0) + '.dys')
# r = [0.3, 0.25, 0.2, 0.15, 0.1, 0.05, 0.04, 0.03, 0.02, 0.01, 0.005, 0.003, 0.002, 0.001]
# r = [0.01]
# nf = np.zeros(len(r))
# for idx, el in enumerate(r):
#     nf[idx] = np.log(hausdorff(idos, el))


# from scipy.stats import linregress
# slope, intercept, r_value, p_value, std_err =linregress(np.log(idos[0]/r/1e9),nf)
# plt.scatter(np.log(idos[0]/r/1e9), nf)
# plt.xlabel('log(f0/r)')
# plt.ylabel('log(N(r)')
# plt.show()
# print(r_value**2, slope)
# print(hausdorff(idos, 0.8))

# plot_full_gap_evolution()
# plt.clf()
# gap_evolution()

# for i in range(100):
#     bragg_sub(i)
#     plt.savefig('sub_' +str(i) +'.png', png=200)
#     plt.close()

# plot_peaks_idos()
# plt.plot(2*np.pi*np.arange(0, len(dys))/lattice_const, dys)
# plt.plot(np.linspace(0.01*2, 2*377, len(struct)), struct[:len(struct)])
# plt.plot(2*np.pi*np.arange(0, len(coef))/lattice_const, coef[:len(coef)])
# plt.plot(k, dispersion(k))
# plt.ylim([10, 40])
# plt.xlim([0,1.75e8])

# plot_peaks_idos()


#
#
# peaks = find_peaks(coef, 0.7)[0]
# print(peaks)
# #
# plt.scatter(x_coef[89], coef[89], color='r')
# plt.vlines(x_coef[peaks], 0, 50)
# plt.scatter(x, dys)
# plt.scatter(x, phason_dys)
# plt.plot(x_coef, coef)
# plt.plot(x, phason_fft)
