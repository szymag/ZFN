import math
import unittest
import  numpy as np
from src.eig_problem.FFTfromFile import FFTfromFile
from src.eig_problem.MacierzDoZagadnienia import MacierzDoZagadnienia
from src.eig_problem.WektorySieciOdwrotnej import WektorySieciOdwrotnej
from src.eig_problem.ParametryMaterialowe import ParametryMaterialowe

class MacierzDoZagadnieniaTestCases(unittest.TestCase):

    def testMatrixtoEignevalue_TheImpact(self):
        tmp = MacierzDoZagadnienia('ff=0.5.txt', 11, 11, np.array([1e-9, 0]), 400e-9, 400e-9,
                                   1.752e6, 0.484e6, 1.09e-17, 5.84e-17,
                                   20e-9, 0, 0.1/(4e-7 * np.pi)).fill_matrix().view(float)
        tmp1 = np.loadtxt('matrix_to_eig_TheImpact_11_11_q=[1e-9,0].txt')
        np.testing.assert_array_almost_equal(tmp, tmp1, decimal=5)

        tmp = MacierzDoZagadnienia('ff=0.5.txt', 3, 3, np.array([0, 1e-9]), 400e-9, 400e-9,
                                   1.752e6, 0.484e6, 1.09e-17, 5.84e-17,
                                   20e-9, 0, 0.1/(4e-7 * np.pi)).fill_matrix().view(float)
        tmp1 = np.loadtxt('matrix_to_eig_TheImpact_3_3_q=[0,1e-9].txt')
        np.testing.assert_array_almost_equal(tmp, tmp1, decimal=5)

        tmp = MacierzDoZagadnienia('ff=0.5.txt', 11, 11, np.array([1e-9, 0]), 400e-9, 400e-9,
                                   1.752e6, 0.484e6, 1.09e-17, 5.84e-17,
                                   20e-9, 0, 0.1/(4e-7 * np.pi)).fill_matrix().view(float)
        tmp1 = np.loadtxt('matrix_to_eig_TheImpact_11_11_q=[1e-9,0].txt')
        np.testing.assert_array_almost_equal(tmp, tmp1, decimal=5)

        tmp = MacierzDoZagadnienia('ff=0.5.txt', 11, 11, np.array([0, 20e-9]), 400e-9, 400e-9,
                                   1.752e6, 0.484e6, 1.09e-17, 5.84e-17,
                                   20e-9, 0, 0.1/(4e-7 * np.pi)).fill_matrix().view(float)
        tmp1 = np.loadtxt('matrix_to_eig_TheImpact_11_11_q=[0,20e-9].txt')
        np.testing.assert_array_almost_equal(tmp, tmp1, decimal=5)

if __name__ == '__main__':
    unittest.main()