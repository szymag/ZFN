# -*- coding: utf-8 -*-

import os
import sys
import unittest

def create_dir(path):
    try:
        os.mkdir(path)
    except OSError:
        print("%s already exists!" % path, file=sys.stderr)


def set_up():
    """Create files and directories for tests (e.g. exported data)"""
    create_dir(os.path.join(os.path.abspath('.'), 'tst', 'tmp'))
    create_dir(os.path.join(os.path.abspath('.'), 'log'))

if __name__ == '__main__':
    set_up()
    print(set_up())
    suite = unittest.TestLoader().discover('.')
    unittest.TextTestRunner(verbosity=2).run(suite)

