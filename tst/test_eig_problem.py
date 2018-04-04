import unittest
import numpy as np
from src.eig_problem.EigenMatrix import EigenMatrix
from src.eig_problem.EigenValueProblem import EigenValueProblem2D

from src.io.DataReader import load_yaml_file


class TheImpactTestCases(unittest.TestCase):
    """Tests for structure in paper "The impact of the lattice symmetry and
       the inclusion shape on the spectrum of 2D magnonic crystals"""
    loaded_data = load_yaml_file("Parameter_for_TheImpact.yaml")
    Fe = loaded_data['material_parameters']['Fe']
    Ni = loaded_data['material_parameters']['Ni']
    dim_sys = loaded_data['system_dimensions']
    # TODO: remove unnecessary parameters
    H0 = 0.1 / (4e-7 * np.pi)
    gamma = 176e9
    mu0H0 = 0.1
    rec_vector_x = 7
    rec_vector_y = 7
    coefficient_from_file = 'ff=0.5.dat'

    tested_case_1 = EigenMatrix(EigenMatrix.ReciprocalVectorGrid(11, 11), np.array([1e-9, 0]),
                                "Parameter_for_TheImpact.yaml", 'Fe', 'Ni')
    tested_case_2 = EigenMatrix(EigenMatrix.ReciprocalVectorGrid(3, 3), np.array([0, 1e-9]),
                                "Parameter_for_TheImpact.yaml", 'Fe', 'Ni')
    tested_case_3 = EigenMatrix(EigenMatrix.ReciprocalVectorGrid(11, 11), np.array([1e-9, 0]),
                                "Parameter_for_TheImpact.yaml", 'Fe', 'Ni')
    tested_case_4 = EigenMatrix(EigenMatrix.ReciprocalVectorGrid(11, 11), np.array([0, 20e-9]),
                                "Parameter_for_TheImpact.yaml", 'Fe', 'Ni')

    tested_case_5 = EigenValueProblem2D('x', "Parameter_for_TheImpact.yaml")

    @staticmethod
    def test_eigen_matrix():
        tested_data_1 = TheImpactTestCases.tested_case_1.generate_and_fill_matrix().view(float)
        tested_data_2 = TheImpactTestCases.tested_case_2.generate_and_fill_matrix().view(float)
        tested_data_3 = TheImpactTestCases.tested_case_3.generate_and_fill_matrix().view(float)
        tested_data_4 = TheImpactTestCases.tested_case_4.generate_and_fill_matrix().view(float)

        loaded_data_to_compare_1 = np.loadtxt('matrix_to_eig_TheImpact_11_11_q=[1e-9,0].dat')
        loaded_data_to_compare_2 = np.loadtxt('matrix_to_eig_TheImpact_3_3_q=[0,1e-9].dat')
        loaded_data_to_compare_3 = np.loadtxt('matrix_to_eig_TheImpact_11_11_q=[1e-9,0].dat')
        loaded_data_to_compare_4 = np.loadtxt('matrix_to_eig_TheImpact_11_11_q=[0,20e-9].dat')

        np.testing.assert_array_almost_equal(loaded_data_to_compare_1, tested_data_1, decimal=5)
        np.testing.assert_array_almost_equal(loaded_data_to_compare_2, tested_data_2, decimal=5)
        np.testing.assert_array_almost_equal(loaded_data_to_compare_3, tested_data_3, decimal=5)
        np.testing.assert_array_almost_equal(loaded_data_to_compare_4, tested_data_4, decimal=5)

    @staticmethod
    def test_eigen_vector():
        tested_case_5 = TheImpactTestCases.tested_case_5.calculate_eigen_vectors().view(float)

        loaded_data_to_compare_5 = np.loadtxt('vec.vec')

        np.testing.assert_array_almost_equal(tested_case_5, loaded_data_to_compare_5)

    @staticmethod
    def test_eig_frequency():
        test_val = np.array(
            [6200586979.692609, 7153839021.040824, 8149434867.91317, 8907414200.202871, 9771172581.694923,
             10249361075.222427, 10275644638.097378, 11961621938.66819, 12102555201.448374, 12218332063.591162,
             12578738960.802431, 12898592836.702763, 12923264289.743351, 13294709487.00567, 13644773768.307749,
             13914853923.230846, 13969207378.576288, 14421398546.441845, 14941243792.562695, 15166512548.983387,
             15781587858.582544, 16552577232.56582, 16769540044.24838, 17204614743.31105, 17269743467.831867,
             17520671539.456398, 17996779666.97687, 18597256883.51367, 19032483290.52158, 19280974273.044727,
             19549741705.78757, 20231896923.40249, 20643618590.448513, 21290731270.40315, 22067497140.468887,
             22425767607.691525, 22640067121.00172, 23585341173.186714, 24229395632.622513, 24480004996.359802,
             25090300474.470634, 25623554039.664497, 25676348194.024654, 26805600771.145847, 27268152868.526615,
             28783565233.363815, 28792615746.10499, 30688563763.639248, 31036486325.37612])

        np.testing.assert_array_almost_equal(EigenValueProblem2D('x',
                                                                 "Parameter_for_TheImpact.yaml").
                                             calculate_eigen_frequency([1e-9, 0]), test_val, decimal=5)


if __name__ == '__main__':
    unittest.main()
