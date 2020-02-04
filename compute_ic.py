#!/usr/bin/python
import os

import matplotlib.pyplot as plt
import numpy as np

from post_process import load


# from scipy.stats import iqr


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
        self.capacity, self.number_of_bins, self.p0_integral = self.calculate_ic()
        self.bins = self.calculate_bins(num_bins=self.number_of_bins)

    def calculate_bins(self, num_bins=100):
        count, self_bins = np.histogram(self.self_output, bins=self.estimator, density=True)
        count, foreign_bins = np.histogram(self.foreign_output, bins=self.estimator, density=True)

        bins = np.linspace(min(self_bins), max(foreign_bins), num=num_bins)

        return bins

    # def kde_plot(self):
    #     sns.distplot(self.foreign_output, hist=True, kde=True,
    #                  bins=self.bins, label="Lf Ls")
    #     sns.distplot(self.self_output, hist=True, kde=True,
    #                  bins=self.bins, label="Ls")
    #     plt.legend()
    # sns.distplot(self.foreign_output, hist=True, kde=True,
    #              bins=self.bins, color='darkblue',
    #              hist_kws={'edgecolor': 'black'},
    #              kde_kws={'linewidth': 4})

    def count_cn(self, bin_locations):
        count_cn, bins = np.histogram(self.foreign_output, bins=bin_locations, density=True)

        return count_cn

    def cn_mean_iqr(self):
        mean = np.mean(self.foreign_output)
        # cn_iqr = iqr(self.foreign_output)

        return mean  # , cn_iqr

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
        # dn_iqr = iqr(self.self_output)
        return mean  # , dn_iqr

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

    # def kullback_leibler(self):
    #     count_cn = self.count_cn()
    #     count_dn = self.count_dn()
    #
    #     bin_width = self.bins[1] - self.bins[0]
    #     kl_cn_dn = np.nan_to_num(np.log2(count_cn/count_dn))
    #     kl_cn_dn = np.trapz(count_cn * np.nan_to_num(np.log2(count_cn/count_dn)), dx=bin_width)
    #     return kl_cn_dn

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
            # print("p(O) integral = " + str(p_0_integral))

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

        return C, number_of_bins, p_0_integral

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

    # def plot_histograms(self):
    #     count_C0, bins_C0, _ = plt.hist(self.foreign_output, bins='auto', align='mid', normed=True,
    #                                     label='P(C{0})'.format(self.num_steps - 1))
    #     print("C bins = " + str(bins_C0))
    #     print("p(CN) integral = " + str(np.trapz(count_C0, dx=bins_C0[1]-bins_C0[0])))
    #
    #     count_D0, bins_D0, _ = plt.hist(self.self_output, bins_C0, align='mid', normed=True,
    #                                     label='P(D{0})'.format(self.num_steps - 1))
    #     print("D bins = " + str(bins_D0))
    #     print("p(DN) integral = " + str(np.trapz(count_D0, dx=bins_D0[1]-bins_D0[0])))
    #
    #     binwidth_C0 = [bins_C0[i+1] - bins_C0[i] for i in range(100)]
    #     print("binwidth = " + str(binwidth_C0))
    #
    #     # Calculating Information capacity
    #
    #     p_O = 0.5 * (count_C0 + count_D0)
    #
    #     print("p(0) integral = " + str(np.trapz(p_O, dx=binwidth_C0[0])))
    #
    #     term_1_C0 = 0.5 * count_C0 * np.nan_to_num(np.log2(count_C0/p_O))
    #     term_2_D0 = 0.5 * count_D0 * np.nan_to_num(np.log2(count_D0/p_O))
    #     C = np.trapz(term_1_C0 + term_2_D0, dx=binwidth_C0[0])
    #     print("C = " + str(C))
    #
    #     np.savetxt("IC", [C], fmt="%f")
    #     np.savetxt("num_bins", [self.num_bins], fmt="%f")
    #     np.savetxt("binwidth", [binwidth_C0[0]], fmt="%f")


def check_binning():
    foreign_output = np.loadtxt("L_self/output")
    foreign_output_end_step = np.loadtxt("3_step_end_step/L_self/output")
    foreign, bins, _ = plt.hist(foreign_output, bins=100, normed=True, label="P(f)")
    foreign_end_step, bins, _ = plt.hist(foreign_output_end_step, bins=100, normed=True, label="P(f) end step")
    plt.legend()

# if __name__ == "__main__":
#     ic = InformationCapacity()
#     ic.calculate_ic()

    # ic.alternate_calculate_ic()
    # check_binning()
    # plt.savefig("output_histograms_test.pdf", format='pdf')
    # ic.plot_histograms()
    # plt.xlim(0,600)
    # plt.legend()
    # plt.savefig("output_histograms.pdf", format='pdf')


