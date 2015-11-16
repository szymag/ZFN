import math
import unittest

from src.eig_problem.MacierzDoZagadnienia import MacierzDoZagadnienia


class Test(unittest.TestCase):
    q = MacierzDoZagadnienia(5)

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

    def test_funkcja_c(self):
        self.assertEqual(self.q.funkcja_c((0, 0), (0, 0), '+'), 1)
        self.assertAlmostEqual(self.q.funkcja_c((1, 1), (0, 0), '+'), 0.9999999972, delta=0.0000000001)
        self.assertAlmostEqual(self.q.funkcja_c((0, 0), (0, 10000000), '-'), 0.98039471957616153, delta=0.0000000000001)
        self.assertAlmostEqual(self.q.funkcja_c((1000, 0), (0, 0), '+'), 0.999998, delta=0.00000001)
        self.assertAlmostEqual(self.q.funkcja_c((2.094395102393195629e+08, 0), (0, 0), '+'), 0.71633974326141991,
                               delta=0.00000001)

    def test_wspolczynnik(self):
        self.assertAlmostEqual(self.q.wspolczynnik((0, 0), (0, 0)), 994611.1809, delta=0.0001)

    def test_czwarte_wyrazenie(self):
        self.assertEqual(self.q.czwarte_wyrazenie((3, 3), (3, 3)), 0)
        self.assertEqual(self.q.czwarte_wyrazenie((3, -3), (3, -3)), 0)


if __name__ == '__main__':
    unittest.main()
