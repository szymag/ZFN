import numpy as np


class FFTzPliku:
    def __init__(self, filepath=None):
        self.tablica = []

    def load_table_from(self, filepath=None):
        if filepath is not None:
            self.tablica = np.loadtxt(filepath).view(complex)
            return self.tablica


    def pr(self):
        print(self.tablica)


if __name__ == "__main__":
    q = FFTzPliku()

    q.pr()
