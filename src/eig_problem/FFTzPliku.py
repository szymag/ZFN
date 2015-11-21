from numpy import genfromtxt


class FFTzPliku:
    def __init__(self):
        self.tablica = genfromtxt('fft1.txt', dtype=complex)

    def pr(self):
        print(self.tablica)


if __name__ == "__main__":
    q = FFTzPliku()

    q.pr()
