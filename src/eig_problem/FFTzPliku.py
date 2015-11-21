from numpy import genfromtxt


class FFTzPliku:
    def __init__(self, filename=None):
        if filename is not None:
            self.tablica = genfromtxt(filename, dtype=complex)
        else:
            self.tablica = []

    def pr(self):
        print(self.tablica)


if __name__ == "__main__":
    q = FFTzPliku()

    q.pr()
