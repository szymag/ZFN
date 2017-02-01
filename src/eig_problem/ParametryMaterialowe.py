from numpy import pi


class ParametryMaterialowe:
# Materials

    # Cobalt 1.445e6 2.287e-17
    MoCo = 1.15e6
    lCo =3.47e-17
    # Permalloy
    MoPy = 0.658e6
    lPy = 4.08e-17
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
    mu0H0 = 0.2
    gamma = 194.6e9
    H0 = mu0H0 / (4e-7 * pi)
    mu0 = 4e-7 * pi
    d = 30e-9 # thickness of material
    a = 500e-9
    x = 0 # position of calculation of dispersion in z-direction
    #r = 14e-9 # radius of inclusion, only for DFT coefficient

# Material parameters
    # inclusion
    MoA = MoCo
    lA = lCo
    # matrix
    MoB = MoPy
    lB = lPy

    def __init__(self):
        pass
