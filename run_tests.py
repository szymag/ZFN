# -*- coding: utf-8 -*-

import unittest
import sys

if __name__ == '__main__':
    suite = unittest.TestLoader().discover('.')
    unittest.TextTestRunner(verbosity=2).run(suite)

