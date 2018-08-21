#!/usr/bin/python
import os

# import matplotlib.pyplot as plt
import numpy as np

from post_process import load


def check_columns(file_path):
    column_names = load(file_path + "column_names")[0].split()
    if "Lf" in column_names[-1] or "Ls" in column_names[-1]:
        self_foreign = True
    else:
        self_foreign = False

    return self_foreign


class PlotOutputHist(object):
    def __init__(self):
        self.column_names = load("sample_0/column_names")[0].split()
        print(self.column_names)

    def compute_output(self):
        output_array = []
        self_foreign = check_columns("sample_0/")
        count = 0
        for i in range(1000):
            try:
                trajectory = np.loadtxt("sample_{0}/mean_traj".format(i))
                if len(trajectory.shape) == 1:
                    if "Ls_Lf" in os.getcwd():
                        if self_foreign:
                            if i == 0:
                                print(self.column_names[-1] + " + " + self.column_names[-2])
                            output_array.append(trajectory[-1] + trajectory[-2])
                        else:
                            if i == 0:
                                print(self.column_names[-1])
                            output_array.append(trajectory[-1])
                    else:
                        if i == 0:
                            # print("Only Self")
                            print(self.column_names[-1])
                        output_array.append(trajectory[-1])
                else:
                    output_array.append(trajectory[-1, -1])
            except:
                print("Error at sample_{0}".format(i))
                count += 1
        print(str(count))
        np.savetxt("output", output_array, fmt='%f')
        return output_array

    # def plot_hist(self):
    #     output = self.compute_output()
    #     count, bins, _ = plt.hist(output, 100, align='mid', normed=True, label='P(O)')
    #     return count


class PlotLigandHist(object):
    def __init__(self):
        self.ligand_concentration = np.loadtxt("Ligand_concentrations")
        self.num_bins = 100

    # def plot_hist(self):
    #     count, bins, _ = plt.hist(self.ligand_concentration, self.num_bins, align='mid', normed=True, label='P(Lf)')

        # mu = 6.0
        # sigma = 1.0
        #
        # x = np.linspace(min(bins), max(bins), 10000)
        # pdf = np.exp(-(np.log(x) - mu) ** 2 / (2 * sigma ** 2)) / (x * sigma * np.sqrt(2 * np.pi))
        # plt.plot(x, pdf, linewidth=2, color='r', label='Log Normal')


if __name__ == "__main__":

    output = PlotOutputHist()
    output.compute_output()
