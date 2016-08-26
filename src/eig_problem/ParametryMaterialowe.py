from numpy import pi

class ParametryMaterialowe:
# Materials
    # Cobalt
    MoCo = 1.445e6
    lCo = 2.287e-17
    # Permalloy
    MoPy = 0.886e6
    lPy = 2.8e-17
    # Iron
    MoFe = 1.752e6
    lFe = 1.09e-17
    # Nickel
    MoNi = 0.484e6
    lNi = 5.84e-17

# System parameters
    mu0H0 = 0.1
    gamma = 176e9
    H0 = mu0H0 / (4e-7 * pi)
    mu0 = 4e-7 * pi
    d = 30e-9 # thickness of material

    x = 0 # position of calculation of dispersion in z-direction
    r = 14e-9 # radius of inclusion, only for DFT coefficient

# Material parameters
    # inclusion
    MoA = MoPy
    lA = lPy
    # matrix
    MoB = MoCo
    lB = lCo

# Program parameters
    # input FFT coefficient

    # amount of reciprocal vectors

    # name of output file

    def __init__(self):
        pass
