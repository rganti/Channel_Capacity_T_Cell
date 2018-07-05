#!/usr/bin/python
import matplotlib.pyplot as plt
import numpy as np


class InformationCapacity(object):

    def __init__(self, num_steps, self_directory="./", foreign_directory="./", estimator='fd', limiting='foreign'):
        self.num_steps = num_steps
        self.foreign_output = np.loadtxt(foreign_directory + "output")
        self.foreign_ligand = np.loadtxt(foreign_directory + "Ligand_concentrations")
        self.self_output = np.loadtxt(self_directory + "output")
        self.self_ligand = np.loadtxt(self_directory + "Ligand_concentrations")

        if limiting == 'foreign':
            self.limiting_output = self.foreign_output
        else:
            self.limiting_output = self.self_output

        self.bins = self.calculate_bins(estimator)

    def calculate_bins(self, estimator):
        count, bins = np.histogram(self.limiting_output, bins=estimator, density=True)
        return bins

    def count_cn(self):
        count_cn, bins = np.histogram(self.foreign_output, bins=self.bins, density=True)
        return count_cn

    def plot_cn(self):
        count, bins, _ = plt.hist(self.foreign_output, bins=self.bins, align='mid', normed=True,
                                  label='P(C{0})'.format(self.num_steps - 1))

    def plot_lf(self, bins):
        count, bins_lf, _ = plt.hist(self.foreign_ligand, bins=bins, align='mid', normed=True,
                                     label='P(Lf)')

    def count_dn(self):
        count_dn, bins = np.histogram(self.self_output, bins=self.bins, density=True)
        return count_dn

    def plot_dn(self):
        count, bins, _ = plt.hist(self.self_output, bins=self.bins, align='mid', normed=True,
                                  label='P(D{0})'.format(self.num_steps - 1))

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


    def calculate_ic(self):
        count_cn = self.count_cn()
        count_dn = self.count_dn()
        bin_width = self.bins[1] - self.bins[0]

        p_O = 0.5 * (count_cn + count_dn)

        print("p(0) integral = " + str(np.trapz(p_O, dx=bin_width)))

        term_1_c0 = 0.5 * count_cn * np.nan_to_num(np.log2(count_cn / p_O))
        term_2_d0 = 0.5 * count_dn * np.nan_to_num(np.log2(count_dn / p_O))
        C = np.trapz(term_1_c0 + term_2_d0, dx=bin_width)
        print("C = " + str(C))
        return C

    def alternate_calculate_ic(self):
        count_cn = self.count_cn()
        count_dn = self.count_dn()
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


if __name__ == "__main__":
    num_steps = 2

    ic = InformationCapacity(num_steps)
    ic.calculate_ic()
    ic.alternate_calculate_ic()
    # check_binning()
    # plt.savefig("output_histograms_test.pdf", format='pdf')
    # ic.plot_histograms()
    # plt.xlim(0,600)
    # plt.legend()
    # plt.savefig("output_histograms.pdf", format='pdf')


