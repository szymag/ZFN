import numpy as np
import multiprocessing as mp
import matplotlib.pyplot as plt

class Chain_Generation:

    def __init__(self, num_blocks, repeat):
        self.num_blocks = num_blocks
        self.repeat = repeat

    def fib_number(self):
        n = self.num_blocks
        i = h = 1
        j = k = 0
        while (n > 0):
            if (n % 2 == 1):
                t = j * h
                j = i * h + j * k + t
                i = i * k + t
            t = h * h
            h = 2 * k * h + t
            k = k * k + t
            n = int(n / 2)
        return j

    def fibonacci_seq(self):
        seq1 = [1]
        seq2 = [0]
        seq = seq2 + seq1
        for i in range(self.num_blocks - 2):
            seq = seq2 + seq1
            seq1 = seq2
            seq2 = seq
        return np.repeat(seq, self.repeat)

    def periodic_seq(self):
        seq = np.zeros(self.fib_number())
        seq[::2] += 1
        return np.repeat(seq, self.repeat)
