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


if __name__ == '__main__':
    unittest.main()
