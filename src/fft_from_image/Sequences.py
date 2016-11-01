from src.fft_from_image.ChainGeneration import ChainGeneration
import numpy as np


class ThueMorse(ChainGeneration):
    def __init__(self, repeat, tm_num):
        ChainGeneration.__init__(self, repeat)
        self.tm_num = tm_num

    @staticmethod
    def tm_construct(seq):
        return [(i + 1) % 2 for i in seq]

    def sequence(self):
        seq = [0]
        for i in range(self.tm_num):
            seq = seq + self.tm_construct(seq)
        return np.repeat(seq, self.repeat)


class Fibonacci(ChainGeneration):
    def __init__(self, repeat, fib_num):
        ChainGeneration.__init__(self, repeat)
        self.fib_num = fib_num

    def fib_number(self):
        n = self.fib_num
        i = h = 1
        j = k = 0
        while n > 0:
            if n % 2 == 1:
                t = j * h
                j = i * h + j * k + t
                i = i * k + t
            t = h * h
            h = 2 * k * h + t
            k = k * k + t
            n = int(n / 2)
        return j

    def sequence(self):
        seq1 = [1]
        seq2 = [0]
        seq = seq2 + seq1
        for i in range(self.fib_num - 2):
            seq = seq2 + seq1
            seq1 = seq2
            seq2 = seq
        return np.repeat(seq, self.repeat)


class Periodic(ChainGeneration):
    def __init__(self, repeat, num):
        ChainGeneration.__init__(self, repeat)
        self.num = num

    def sequence(self):
        seq = np.zeros(self.num)
        seq[::2] += 1
        #seq %= 2
        return np.repeat(seq, self.repeat)

class Heated(ChainGeneration):
    def __init__(self, repeat):
        ChainGeneration.__init__(self, repeat)

    def cos_sequence(self):
        return (np.cos(np.linspace(0, 2 * np.pi, self.repeat)) + 1) / 2