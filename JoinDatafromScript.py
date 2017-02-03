import numpy as np

class JoinFilefromScript:
    def __init__(self, begin_num, end_num, first_part_of_name):
        self.begin_num = begin_num
        self.end_num = end_num
        self.first_part_of_name = first_part_of_name


    def load_data_and_join(self):
        joined_data = []
        for i in range(self.begin_num, self.end_num):
            joined_data.append(np.loadtxt(self.first_part_of_name + str(i) + ".dat"))
            np.savetxt('cni0.5.dat', np.delete(np.transpose(np.array(joined_data)), 0, 0))


class FieldValue:
    def __init__(self, loaded_file, freq_cut):
        self.freq_vs_field = np.loadtxt(loaded_file)
        self.freq_cut = freq_cut

    def place_nearest_freq(self, branch_modes):
        return np.argmin(abs(branch_modes - self.freq_cut))

    def tmp(self, num_modes):
        tmp = []
        for i in range(num_modes + 1):
            tmp.append(self.place_nearest_freq(self.freq_vs_field[i]))
        np.savetxt('fcni0.5.dat', tmp)
        return tmp


#q = JoinFilefromScript(1, 400, '0.5cdys_')
#q.load_data_and_join()
q = FieldValue('cni0.5.dat', 5e9).tmp(6)
print(q)