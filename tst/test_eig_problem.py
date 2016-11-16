import math
import unittest
import  numpy as np
from src.eig_problem.FFTfromFile import FFTfromFile
from src.eig_problem.MacierzDoZagadnienia import MacierzDoZagadnienia
from src.eig_problem.WektorySieciOdwrotnej import WektorySieciOdwrotnej

class MacierzDoZagadnieniaTestCases(unittest.TestCase):
    parameters = np.genfromtxt('set_of_parameters.dat', dtype=None, delimiter=', ')
    #tmp = MacierzDoZagadnienia(parameters)
    print(parameters[0][0].astype(str))

    def build_objects(self):
        pass
        #for i in enumerate(parameters)

    def test_pole_wymiany_II(self):
        ans_pole_wymiany_II = np.genfromtxt('answers_to_pole_wymiany_II.dat', dtype=None, delimiter=' ')
        pass

if __name__ == '__main__':
    unittest.main()