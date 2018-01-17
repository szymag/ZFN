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
    MoNi = 0.484e6 # saturation magnetization
    lNi = 5.84e-17 # exchange constant
    # YIG
    MoY = 0.194e6
    lY = 1.7e-16
    # CoFeB
    MoCoFeB = 1.25e6
    lCoFeB = 1.53e-17
# System parameters
    mu0H0 = None
    gamma = 176e9
    mu0 = 4e-7 * pi
    H0 = None
    d = 60e-9 # thickness of material
    a = 1100e-9
    x = 0 # position of calculation of dispersion in z-direction
    angle = None
    #r = 14e-9 # radius of inclusion, only for DFT coefficient

# Material parameters
    # inclusion
    MoA = 0.915*MoNi
    lA = 0.915*lNi
    # matrix
    MoB = 0.909*MoNi
    lB = 0.909*lNi

    def __init__(self):
        pass
