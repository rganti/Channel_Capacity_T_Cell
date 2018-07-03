#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt


class PlotOutputHist(object):
    def __init__(self):
        pass

    def compute_output(self):
        output_array = []
        for i in range(1000):
            try:
                trajectory = np.loadtxt("sample_{0}/mean_traj".format(i))
                if len(trajectory.shape) == 1:
                    output_array.append(trajectory[-1])
                else:
                    output_array.append(trajectory[-1, -1])
            except:
                print("Error at sample_{0}".format(i))
        np.savetxt("output", output_array, fmt='%f')
        return output_array

    def plot_hist(self):
        output = self.compute_output()
        count, bins, _ = plt.hist(output, 100, align='mid', normed=True, label='P(O)')
        return count


class PlotLigandHist(object):
    def __init__(self):
        self.ligand_concentration = np.loadtxt("Ligand_concentrations")
        self.num_bins = 100

    def plot_hist(self):
        count, bins, _ = plt.hist(self.ligand_concentration, self.num_bins, align='mid', normed=True, label='P(Lf)')

        # mu = 6.0
        # sigma = 1.0
        #
        # x = np.linspace(min(bins), max(bins), 10000)
        # pdf = np.exp(-(np.log(x) - mu) ** 2 / (2 * sigma ** 2)) / (x * sigma * np.sqrt(2 * np.pi))
        # plt.plot(x, pdf, linewidth=2, color='r', label='Log Normal')


if __name__ == "__main__":

    output = PlotOutputHist()
    output.plot_hist()
    plt.savefig("output.pdf", format='pdf')

    ligand = PlotLigandHist()
    ligand.plot_hist()
    plt.legend()

    plt.savefig("lf_output.pdf", format='pdf')
    plt.close()