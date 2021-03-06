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

    def sequence_generator(self):
        seq1 = [1]
        seq2 = [0]
        seq = seq2 + seq1
        for i in range(self.fib_num - 2):
            seq = seq2 + seq1
            seq1 = seq2
            seq2 = seq
        return np.array(seq)

    def sequence(self):
        return np.repeat(self.sequence_generator(), self.repeat)


class Periodic(ChainGeneration):
    def __init__(self, repeat, num):
        ChainGeneration.__init__(self, repeat)
        self.num = num

    def sequence_generator(self):
        seq = np.zeros(self.num)
        seq[::2] += 1
        return seq

    def sequence(self):
        return np.repeat(self.sequence_generator(), self.repeat)


class Random(ChainGeneration):
    def __init__(self, repeat, num, stripes1_count):
        ChainGeneration.__init__(self, repeat)
        self.num = num
        self.stripes1_count = stripes1_count

    def sequence_generator(self):
        seq = np.zeros(self.num)
        seq[:self.stripes1_count] += 1
        return np.random.permutation(seq)

    def sequence(self):
        return np.repeat(self.sequence_generator(), self.repeat)


class Heated(ChainGeneration):
    def __init__(self, repeat):
        ChainGeneration.__init__(self, repeat)

    def cos_sequence(self):
        return (np.cos(np.linspace(0, 2 * np.pi, self.repeat)) + 1) / 2


class Custom(ChainGeneration):
    def __init__(self, file_name, repeat=1):
        ChainGeneration.__init__(self, repeat)
        self.tmp = np.transpose(np.loadtxt(file_name))[-1]
        self.data = (self.tmp - np.min(self.tmp))
        self.data /= np.max(self.data)

    def sequence(self):
        return self.data


class Phason:
    def __init__(self, sequence_type, repeat, num, phason_parameter):
        if sequence_type == 'F':
            self.f = Fibonacci(1, num)
            self.seq = self.f.sequence()
        elif sequence_type == 'P':
            self.p = Periodic(1, num)
            self.seq = self.p.sequence()
        elif sequence_type == 'R':
            struct = Fibonacci(1, num).fib_number()
            stripes1 = Fibonacci(1, num - 2).fib_number()
            self.p = Random(1, struct, stripes1) # randomized Fibonacci
            self.seq = self.p.sequence()
        else:
            raise ValueError('No more types supported at the moment')
        self.repeat = repeat
        self.len = len(self.seq)
        self.where_one = self.find_all_phasons(self.seq)
        self.phason_parameter = phason_parameter
        self.sequence_type  = sequence_type
        if phason_parameter < 1:
            self.phasons_count = int(phason_parameter * len(self.where_one))
        else:
            self.phasons_count = phason_parameter

    def find_all_phasons(self, seq):
        a = np.argwhere(seq == 1).T[0]
        b = np.concatenate((np.diff(a), np.array([(self.len - a[-1] + a[0])])))
        return np.compress(np.where(b == 1, 0, 1) == 1, a)

    def sequence_shuffling(self, seq):
        if self.sequence_type == "R":
            phason_pos =  np.argwhere(self.p.sequence_generator() == 1)
            print(phason_pos)
            return self.seq, phason_pos
        else:
            if self.phason_parameter < 1:
                phasons_pos = np.random.permutation(self.find_all_phasons(seq))[0:self.phasons_count]
                seq = self.make_shufling(phasons_pos)
            else:
                collect_phasons = np.zeros(self.phasons_count)
                for i in range(self.phasons_count):
                    phasons_pos = np.random.permutation(self.find_all_phasons(seq))[0]
                    seq = self.make_shufling(phasons_pos)
                    collect_phasons[i] = phasons_pos
                phasons_pos = collect_phasons
            return seq, phasons_pos

    def make_shufling(self, stripe_position):
        seq = self.seq
        seq[(stripe_position + 1) % len(seq)] = 1
        seq[stripe_position] = 0
        return seq

    def sequence(self, seq):
        return np.repeat(seq, self.repeat)


if __name__ == "__main__":
    pass
