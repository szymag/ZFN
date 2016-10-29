from numpy import pi


class ParametryMaterialowe:
    # Materials
    # Cobalt
    MoCo = 1.15e6
    lCo = 3.46e-17
    # Permalloy
    MoPy = 0.658e6
    lPy = 4.05e-17
    # Iron
    MoFe = 1.752e6
    lFe = 1.09e-17
    # Nickel
    MoNi = 0.484e6
    lNi = 5.84e-17
    # YIG
    MoY = 0.194e6
    lY = 1.7e-16

    # System parameters
    mu0H0 = 0.1
    gamma = 176e9
    H0 = mu0H0 / (4e-7 * pi)
    mu0 = 4e-7 * pi
    d = 20e-9  # thickness of material
    a = 400e-9
    x = 0  # position of calculation of dispersion in z-direction

    # Material parameters
    # inclusion
    MoA = MoFe
    lA = lFe
    # matrix
    MoB = MoNi
    lB = lNi

    # Program parameters
    # input FFT coefficient

    # amount of reciprocal vectors

    # name of output file

    def __init__(self):
        pass
