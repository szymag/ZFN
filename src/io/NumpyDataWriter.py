from src.io.DataWriter import DataWriter
import numpy as np


class NumpyDataWriter(DataWriter):
    def __init__(self, pathname):
        super().__init__()
        self.pathname = pathname

    def write(self, data):
        np.savetxt(self.pathname, data)