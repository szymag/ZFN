import math
import unittest

from src.eig_problem.FFTfromFile import FFTfromFile
from src.eig_problem.MacierzDoZagadnienia import MacierzDoZagadnienia
from src.eig_problem.WektorySieciOdwrotnej import WektorySieciOdwrotnej


class TestEigProblem(unittest.TestCase):
    q = MacierzDoZagadnienia(9, 'FFT', 'I')
    w = FFTfromFile(9, 'I')
    e = WektorySieciOdwrotnej(30e-9, 30e-9, 9)
    ee = WektorySieciOdwrotnej(300e-9, 300e-9, 9)

    def test_suma_roznica_wektorow(self):
        self.assertEqual(self.q.suma_roznica_wektorow((0, 0), (0, 0), '+'), (0, 0))
        self.assertEqual(self.q.suma_roznica_wektorow((0, 0), (0, 0), '-'), (0, 0))
        self.assertEqual(self.q.suma_roznica_wektorow((-3, 3), (0, 0), '+'), (-3, 3))
        self.assertEqual(self.q.suma_roznica_wektorow((0, 15), (5, 0), '-'), (-5, 15))
        self.assertEqual(self.q.suma_roznica_wektorow((0, 0), (0, 10), '+'), (0, 10))
        self.assertEqual(self.q.suma_roznica_wektorow((0, 3), (0, 3), '-'), (0, 0))
        self.assertEqual(self.q.suma_roznica_wektorow((4.12, 2), (0, 3120), '+'), (4.12, 3122))
        self.assertEqual(self.q.suma_roznica_wektorow((0, 0), (121340, 0), '+'), (121340, 0))

    def test_norma_wektorow(self):
        self.assertEqual(self.q.norma_wektorow((0, 0), (0, 0), '+'), 0)
        self.assertEqual(self.q.norma_wektorow((1, 1), (0, 0), '-'), math.sqrt(2))
        self.assertEqual(self.q.norma_wektorow((1, 1), (0, 0), '+'), math.sqrt(2))
        self.assertEqual(self.q.norma_wektorow((1, 1), (-1, -1), '-'), math.sqrt(8))
        self.assertEqual(self.q.norma_wektorow((1, 1), (1, -1), '+'), math.sqrt(4))
        self.assertEqual(self.q.norma_wektorow((1, 1), (1, -1), '-'), math.sqrt(4))
        self.assertEqual(self.q.norma_wektorow((100, 11.5), (10, -20), '+'), math.sqrt(110 ** 2 + (-8.5) ** 2))
        self.assertEqual(self.q.norma_wektorow((10, 10), (-50, -50), '+'), math.sqrt(40 ** 2 + 40 ** 2))
        self.assertEqual(self.q.norma_wektorow((1000, 0), (0, 0), '+'), 1000)
