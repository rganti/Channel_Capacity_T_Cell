import numpy as np
from matplotlib import pyplot as plt

from post_process import load


class Check(object):
    def __init__(self, file_path, n):
        self.ordered_network = load(file_path + "ordered_network")
        self.reactions = load(file_path + "kp_steps_{0}.rxn".format(n))
        self.qsub = load(file_path + "qsub.sh")


class CheckEquality(object):
    def __init__(self, file_path_1, file_path_2):
        self.file_1 = load(file_path_1)
        self.file_2 = load(file_path_2)

    def check(self):
        i = 0
        for i in range(len(self.file_1)):
            x = self.file_1[i] == self.file_2[i]

            if not x:
                print("File 1: " + self.file_1[i])
                print("File 2: " + self.file_2[i])

            i += 1

        return x


def plot_eta_previous(file_path, eta_class):
    eta = np.loadtxt(file_path + "eta")

    foreign_column = load(file_path + "Lf/column_names")
    foreign_column_names = foreign_column[0].split()

    self_column = load(file_path + "Ls/column_names")
    self_column_names = self_column[0].split()

    plt.plot(eta_class.rates, np.zeros_like(eta_class.rates) + eta,
             label="{0}/{1}".format(self_column_names[-1], foreign_column_names[-1]),
             linestyle='-', marker='o')


class PlotEta(object):
    def __init__(self, file_path="./"):
        self.eta = np.loadtxt(file_path + "eta")
        self.rates = np.loadtxt(file_path + "rates")

        self.foreign_column = load(file_path + "parameter_0/Lf/column_names")
        self.foreign_column_names = self.foreign_column[0].split()

        self.self_column = load(file_path + "parameter_0/Ls/column_names")
        self.self_column_names = self.self_column[0].split()

    def plot_eta(self):
        plt.plot(self.rates, self.eta,
                 label="{0}/{1}".format(self.self_column_names[-1], self.foreign_column_names[-1]),
                 linestyle='-', marker='o')
