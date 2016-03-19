from numpy import pi


class ParametryMaterialowe:
# Materials
    # Cobalt
    MoCo = 1.145e6
    lCo = 2.3e-17
    # Permalloy
    MoPy = 0.658e6
    lPy = 4.044e-17
    # Iron
    MoFe = 1.752e6
    lFe = 1.09e-17
    # Nickel
    MoNi = 0.484e6
    lNi = 5.84e-17

# System parameters
    mu0H0 = 0.2
    gamma = 194.6e9
    H0 = mu0H0 / (4e-7 * pi)
    mu0 = 4e-7 * pi
    d = 4e-9 # thickness of material
    a = 30e-9 # size of lattice in x-direction
    b = 30e-9 # size of lattice in y-direction
    x = 0 # position of calculation of dispersion in z-direction
    r = 14e-9 # radius of inclusion, only for DFT coefficient

# Material parameters
    # inclusion
    MoA = MoCo
    lA = lCo
    # matrix
    MoB = MoPy
    lB = lPy

# Program parameters
    # input FFT coefficient
    input_fft = 'F7.txt'
    # amount of reciprocal vectors
    ilosc_wektorow = 225
    # name of output file
    outpu_file = 'test.txt'

    def __init__(self):
        pass
