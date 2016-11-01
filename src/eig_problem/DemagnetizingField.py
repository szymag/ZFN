import numpy as np
from src.eig_problem.ParametryMaterialowe import ParametryMaterialowe
from src.eig_problem.WektorySieciOdwrotnej import WektorySieciOdwrotnej
from src.eig_problem.FFTfromFile1D import FFTfromFile1D
from src.eig_problem.ZagadnienieWlasne import ZagadnienieWlasne

class DemagnetizingField:
    def __init__(self, input_fft):
        self.tmp = FFTfromFile1D(input_fft)
        self.ilosc_wektorow = self.tmp.ilosc_wektorow
        self.input_fft = input_fft

    def demagnetizing_field(self):
        a = ZagadnienieWlasne(1, self.input_fft, 'r.xtx').wektory_wlasne()
        print(a)
        return 0

    def inverse_fourier_transform(self):
        pass

if __name__ == "__main__":
    q = DemagnetizingField('p_coef_10*2.txt').demagnetizing_field()