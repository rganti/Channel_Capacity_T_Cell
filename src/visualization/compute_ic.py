import os

import matplotlib.pyplot as plt
import numpy as np

from src.general.directory_handling import load


class InformationCapacity(object):

    def __init__(self, foreign_directory="./", self_directory="./", estimator='fd', limiting='foreign'):
        self.num_steps = 1
        self.foreign_directory = foreign_directory
        self.self_directory = self_directory

        self.foreign_output = np.loadtxt(foreign_directory + "output")
        self.foreign_ligand = np.loadtxt(foreign_directory + "Ligand_concentrations")
        self.self_output = np.loadtxt(self_directory + "output")
        self.self_ligand = np.loadtxt(self_directory + "Ligand_concentrations")

        if os.path.exists(foreign_directory + "sample_0/column_names"):
            print("Loaded foreign column names")
            self.foreign_column = load(foreign_directory + "sample_0/column_names")
            self.foreign_column_names = self.foreign_column[0].split()
        elif os.path.exists(foreign_directory + "column_names"):
            print("Loaded foreign column names")
            self.foreign_column = load(foreign_directory + "column_names")
            self.foreign_column_names = self.foreign_column[0].split()

        if os.path.exists(self_directory + "sample_0/column_names"):
            self.self_column = load(self_directory + "sample_0/column_names")
            self.self_column_names = self.self_column[0].split()
        elif os.path.exists(self_directory + "column_names"):
            self.self_column = load(self_directory + "column_names")
            self.self_column_names = self.self_column[0].split()

        self.estimator = estimator
        self.capacity = self.calculate_ic()

    def calculate_bins(self, num_bins=100):
        count, self_bins = np.histogram(self.self_output, bins=self.estimator, density=True)
        count, foreign_bins = np.histogram(self.foreign_output, bins=self.estimator, density=True)

        bins = np.linspace(min(self_bins), max(foreign_bins), num=num_bins)

        return bins

    def count_cn(self, bin_locations):
        count_cn, bins = np.histogram(self.foreign_output, bins=bin_locations, density=True)

        return count_cn

    def cn_mean_iqr(self):
        mean = np.mean(self.foreign_output)

        return mean

    def plot_cn(self, plot_label=None):
        if os.path.exists(self.foreign_directory + "sample_0/column_names") or \
                os.path.exists(self.foreign_directory + "column_names"):
            if len(self.foreign_column_names) > 1:
                if ("Lf" in self.foreign_column_names[-2] or "Lf" in self.foreign_column_names[-1]) and \
                        ("Ls" in self.foreign_column_names[-2] or "Ls" in self.foreign_column_names[-1]):
                    label = 'P(' + self.foreign_column_names[-2] + " + " + self.foreign_column_names[-1] + ')'
                else:
                    label = 'P(' + self.foreign_column_names[-1] + ')'
            else:
                label = 'P(' + self.foreign_column_names[-1] + ')'
        else:
            label = 'P(C{0} + D{0})'.format(self.num_steps - 1)

        if plot_label:
            label = plot_label

        count, bins, _ = plt.hist(self.foreign_output, bins=self.bins, align='mid', normed=True,
                                  label=label)

    def plot_lf(self, bins):
        count, bins_lf, _ = plt.hist(self.foreign_ligand, bins=bins, align='mid', normed=True,
                                     label='P(Lf)')

    def count_dn(self, bin_locations):
        count_dn, bins = np.histogram(self.self_output, bins=bin_locations, density=True)
        return count_dn

    def dn_mean_iqr(self):
        mean = np.mean(self.self_output)

        return mean

    def plot_dn(self, plot_label=None):
        if os.path.exists(self.self_directory + "sample_0/column_names") or \
                os.path.exists(self.self_directory + "column_names"):
            label = 'P(' + self.self_column_names[-1] + ')'
        else:
            label = 'P(D{0})'.format(self.num_steps - 1)

        if plot_label:
            label = plot_label

        count, bins, _ = plt.hist(self.self_output, bins=self.bins, align='mid', normed=True,
                                  label=label)

    def plot_ls(self, bins):
        count, bins_ls, _ = plt.hist(self.self_ligand, bins=bins, align='mid', normed=True,
                                     label='P(Ls)')

    def compute_bins(self):
        number_of_bins = 50
        bins = 0
        p_0_integral = 0
        while p_0_integral < 0.99:
            bins = self.calculate_bins(num_bins=number_of_bins)
            count_cn = self.count_cn(bins)
            count_dn = self.count_dn(bins)
            bin_width = bins[1] - bins[0]

            p_O = 0.5 * (count_cn + count_dn)

            p_0_integral = np.trapz(p_O, dx=bin_width)

            number_of_bins += 50

        return bins

    def calculate_ic(self):
        number_of_bins = 50
        C = 0
        p_0_integral = 0
        while p_0_integral < 0.99:
            bins = self.calculate_bins(num_bins=number_of_bins)
            count_cn = self.count_cn(bins)
            count_dn = self.count_dn(bins)
            bin_width = bins[1] - bins[0]

            p_O = 0.5 * (count_cn + count_dn)

            p_0_integral = np.trapz(p_O, dx=bin_width)
            print("p(O) integral = " + str(p_0_integral))

            term_1_c0 = 0.5 * count_cn * np.nan_to_num(np.log2(count_cn / p_O))
            term_2_d0 = 0.5 * count_dn * np.nan_to_num(np.log2(count_dn / p_O))
            C = np.trapz(term_1_c0 + term_2_d0, dx=bin_width)
            print("C = " + str(C))

            if p_0_integral == C:
                print("C == P(O) integral: " + str(p_0_integral == C))
                C = 1.00
                print("New C " + str(C))
                break

            number_of_bins += 50

        return C

    def alternate_calculate_ic(self):
        bins = self.calculate_bins(num_bins=500)
        count_cn = self.count_cn(bins)
        count_dn = self.count_dn(bins)
        bin_width = self.bins[1] - self.bins[0]

        p_O = 0.5 * (count_cn + count_dn)

        h_o = -p_O * np.nan_to_num(np.log2(p_O))
        h_o_i_term_1 = 0.5 * count_cn * np.nan_to_num(np.log2(count_cn))
        h_o_i_term_2 = 0.5 * count_dn * np.nan_to_num(np.log2(count_dn))

        C = np.trapz(h_o + h_o_i_term_1 + h_o_i_term_2, dx=bin_width)
        print("C2 = " + str(C))
        return C


def check_binning():
    foreign_output = np.loadtxt("L_self/output")
    foreign_output_end_step = np.loadtxt("3_step_end_step/L_self/output")
    foreign, bins, _ = plt.hist(foreign_output, bins=100, normed=True, label="P(f)")
    foreign_end_step, bins, _ = plt.hist(foreign_output_end_step, bins=100, normed=True, label="P(f) end step")
    plt.legend()



