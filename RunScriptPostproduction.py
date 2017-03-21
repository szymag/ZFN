import numpy as np


class RunScriptPostproduction:
    def __init__(self, begin_num, end_num, first_part_of_name, freq_cut=None):
        self.begin_num = begin_num
        self.end_num = end_num
        self.first_part_of_name = first_part_of_name
        self.freq_cut = freq_cut

    def field_for_given_freq(self, num_modes):
        tmp = []
        for i in range(num_modes + 1):
            tmp.append(self.place_nearest_desirable_freq(self.load_data_and_join()[i]))
        np.savetxt('fcni0.5.dat', tmp)
        return tmp

    def place_nearest_desirable_freq(self, branch_modes):
        return np.argmin(abs(branch_modes - self.freq_cut))

    def load_data_and_join(self):
        joined_data = []
        for i in range(self.begin_num, self.end_num):
            joined_data.append(np.loadtxt(self.first_part_of_name + str(i) + ".dat"))
            np.savetxt('0_fmr.dat', np.delete(np.transpose(np.array(joined_data)), 0, 0))
            #return np.delete(np.transpose(np.array(joined_data)), 0, 0)


if __name__ == "__main__":
    RunScriptPostproduction(1, 400, '0_',).load_data_and_join()
