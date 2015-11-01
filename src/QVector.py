__author__ = 'szymag'

from numpy import cosh, linalg


class Qvector:
    def __init__(self, vector):
        self.vector = vector

    def q(self):
        return cosh(linalg.norm(self.vector))
