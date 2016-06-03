import numpy as np
import multiprocessing as mp


class Chain_Generation:

    def __init__(self, num_blocks, fib_num):
        self.num_blocks = num_blocks
        self.fib_num = fib_num

    def fib_matrix(self, n):
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

    def Fibonacci(self):
        pool = mp.Pool()
        assert len != 0 or len != 1, 'The sequence is too short'
        seq = pool.map(self.fib_matrix, np.arange(self.num_blocks))
        return seq

    def fibonacci_seq(self):
        seq = np.zeros(self.fib_matrix(self.fib_num))
        seq[1] = 1
        for i in list(range(2,self.fib_num)):
            pass


