import unittest
from src.eig_problem.EigenMatrix import EigenMatrix
import numpy as np

class ReciprocalVectorGridTestCases(unittest.TestCase):

    def test_raise_on_invalid_arguments(self):
        with self.assertRaises(ValueError):
            grid = EigenMatrix.ReciprocalVectorGrid(0, 0)
        with self.assertRaises(ValueError):
            grid = EigenMatrix.ReciprocalVectorGrid(-1,1)

    def test_return_correct_vectors_count(self):
        self.assertEqual(EigenMatrix.ReciprocalVectorGrid(5, 5).vectors_count(), 25)
        self.assertEqual(EigenMatrix.ReciprocalVectorGrid(5, 4).vectors_count(), 20)
        self.assertEqual(EigenMatrix.ReciprocalVectorGrid(4, 5).vectors_count(), 20)

    def test_return_correct_value_grid_size(self):
        grid_size = EigenMatrix.ReciprocalVectorGrid(9, 9).coefficient_grid_size()
        self.assertTrue(grid_size[0] % 2 == 1)
        self.assertTrue(grid_size[1] % 2 == 1)
        self.assertEqual(grid_size, (17, 17))

    def test_return_correct_value_sift_to_middle_array(self):
        shift = EigenMatrix.ReciprocalVectorGrid(11, 11).shift_to_middle_of_coeff_array()
        np.testing.assert_array_equal(shift, np.array([10, 10]))
        shift2 = EigenMatrix.ReciprocalVectorGrid(25, 25).shift_to_middle_of_coeff_array()
        np.testing.assert_array_equal(shift2, np.array([24, 24]))