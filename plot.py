#!/usr/bin/python
import argparse

import numpy as np
from matplotlib import pyplot as plt

from two_species import SelfWithForeign


def find_index(dummy_file, string):
    index_string = [string in x for x in dummy_file].index(True)
    return index_string


def load_columns(filename):
    column_names = open(filename, "r").readline().split()
    return column_names


class PlotTwoSpecies(object):

    def __init__(self, file):
        self.data = np.loadtxt(file)
        self.alpha = 0.1
        self.beta = 0.2
        self.A_0 = 800
        self.B_0 = 200

    def plot_concentrations(self):
        plt.plot(self.data[:, 0], self.data[:, 1]/(self.A_0 + self.B_0), label="A", linestyle='-', marker='o')
        plt.plot(self.data[:, 0], self.data[:, 2]/(self.A_0 + self.B_0), label="B", linestyle='-', marker='o')

    def plot_analytical(self):
        x = np.linspace(0,50,500)
        A = (self.A_0 * ((self.beta/(self.alpha + self.beta)) + (self.alpha/(self.alpha + self.beta)) * np.exp(-(self.alpha + self.beta) * x)) +\
            self.B_0 * ((self.beta/(self.alpha + self.beta)) - (self.beta/(self.alpha + self.beta)) * np.exp(-(self.alpha + self.beta) * x)))/(self.A_0 + self.B_0)
        plt.plot(x, A, label="A analytical", linestyle='--')

        B = (self.A_0 * ((self.alpha / (self.alpha + self.beta)) - (self.alpha / (self.alpha + self.beta)) * np.exp(-(self.alpha + self.beta) * x)) +\
            self.B_0 * ((self.alpha / (self.alpha + self.beta)) + (self.beta/ (self.alpha + self.beta)) * np.exp(-(self.alpha + self.beta) * x)))/(self.A_0 + self.B_0)
        plt.plot(x, B, label="B analytical", linestyle='--')

    def save_plot(self):
        self.plot_concentrations()
        self.plot_analytical()
        plt.legend()
        plt.savefig("concentrations.pdf", format='pdf')


class PlotKP(object):

    def __init__(self, filename, num_kp_steps, directory="./"):
        self.column_names = load_columns(directory + "column_names")
        self.num_kp_steps = num_kp_steps

        self.data = np.loadtxt(filename)

        self.time_index = find_index(self.column_names, "time")
        self.time = self.data[:, self.time_index]

        self.Lf_index = find_index(self.column_names, "Lf")
        self.Lf = self.data[:, self.Lf_index]

        self.CN_index = find_index(self.column_names, "C" + str(self.num_kp_steps - 1))
        self.CN = self.data[:, self.CN_index]

        self.Ls_index = find_index(self.column_names, "Ls")
        self.Ls = self.data[:, self.Ls_index]

        self.DN_index = find_index(self.column_names, "D" + str(self.num_kp_steps - 1))
        self.DN = self.data[:, self.DN_index]

    def plot_concentrations(self):
        for i in range(1, len(self.column_names)):
            plt.plot(self.time, self.data[:, i],
                     label=self.column_names[i], linestyle='-', marker='o')

    def error_fraction(self):
        plt.plot(self.time, self.DN/self.CN, label="error rate {0}".format(self.num_kp_steps), linestyle='-', marker='o')
        plt.ylim(0, 0.05)
        plt.legend()
        np.savetxt("error", self.DN/self.CN)
        plt.savefig("error.pdf", format='pdf')

    # def error_fraction_three_steps(self):
    #     if self.num_kp_steps == 3:
    #         C1_index = find_index(self.column_names, "C1")
    #         C1 = self.data[:, C1_index]
    #
    #         D1_index = find_index(self.column_names, "D1")
    #         D1 = self.data[:, D1_index]
    #
    #         C2_index = find_index(self.column_names, "C2")
    #         C2 = self.data[:, C2_index]
    #
    #         D2_index = find_index(self.column_names, "D2")
    #         D2 = self.data[:, D2_index]
    #
    #         plt.plot(self.time, D1/C1, label="D1/C1", linestyle='-', marker='o')
    #         plt.plot(self.time, (D1/C1)**2, label='(D1/C1)^2', linestyle='--')
    #         plt.plot(self.time, D2/C2, label="(D2/C2)", linestyle='-', marker='o')
    #
    #         plt.ylim(0, 0.2)
    #         plt.legend()
    #         plt.savefig("error_2_step.pdf", format='pdf')

    # def plot_analytical(self):
    #     x = np.linspace(0, 50, 500)
    #     A = (self.A_0 * ((self.beta/(self.alpha + self.beta)) + (self.alpha/(self.alpha + self.beta)) * np.exp(-(self.alpha + self.beta) * x)) +\
    #         self.B_0 * ((self.beta/(self.alpha + self.beta)) - (self.beta/(self.alpha + self.beta)) * np.exp(-(self.alpha + self.beta) * x)))/(self.A_0 + self.B_0)
    #     plt.plot(x, A, label="Lf analytical", linestyle='--')
    #
    #     B = (self.A_0 * ((self.alpha / (self.alpha + self.beta)) - (self.alpha / (self.alpha + self.beta)) * np.exp(-(self.alpha + self.beta) * x)) +\
    #         self.B_0 * ((self.alpha / (self.alpha + self.beta)) + (self.beta/ (self.alpha + self.beta)) * np.exp(-(self.alpha + self.beta) * x)))/(self.A_0 + self.B_0)
    #     plt.plot(x, B, label="CO analytical", linestyle='--')

    def save_plot(self):
        self.plot_concentrations()
        # self.plot_analytical()
        # plt.ylim(0, 1)
        plt.legend()
        plt.savefig("concentrations.pdf", format='pdf')
        plt.close()


class PlotSecondOrder(object):

    def __init__(self, filename, directory="./"):
        self.column_names = load_columns(directory + "column_names")
        self.data = np.loadtxt(filename)

        self.time_index = find_index(self.column_names, "time")
        self.time = self.data[:, self.time_index]

        # self.Lf_index = find_index(self.column_names, "Ls")
        # self.Lf = self.data[:, self.Lf_index]
        #
        # self.CN_index = find_index(self.column_names, "D0")
        # self.CN = self.data[:, self.CN_index]

    def plot_concentrations(self):
        # for i in range(1, len(self.column_names)):
        plt.plot(self.time, self.data[:, -1], label=self.column_names[-1], linestyle='-', marker='o')
        plt.plot(self.time, self.data[:, -2], label=self.column_names[-2], linestyle='-', marker='o')

    def plot_analytical(self):
        Lf = self.second_order.n_initial["Lf"]
        R = self.second_order.n_initial["R"]
        kappa = self.second_order.rate_constants.kappa
        k_off = self.second_order.rate_constants.k_off_foreign

        c_o = 0.5 * (-np.sqrt((-Lf - R - (k_off/kappa))**(2) - 4*Lf*R) + Lf + R + (k_off/kappa))
        print("c_o = " + str(c_o))
        plt.plot(self.time, np.zeros_like(self.time) + c_o, label="analytical C0", linestyle='-')

    def save_plot(self):
        self.plot_concentrations()
        plt.legend()
        plt.savefig("concentrations.pdf", format='pdf')
        plt.close()


class MutualInformation(object):
    def __init__(self, kp_data, save=False):
        self.KP_data = kp_data
        self.save = save
        self.KP = SelfWithForeign()

        self.Lf_initial = self.KP.n_initial["Lf"]
        print("Lf_initial = " + str(self.Lf_initial))
        self.Ls_initial = self.KP.n_initial["Ls"]
        print("Ls_initial = " + str(self.Ls_initial))

        # Prior Distributions
        self.p_f = 0.5
        self.p_s = 1 - self.p_f

        # Conditional Distributions

        # self.p_1_given_f = (np.zeros_like(self.KP_data.CN) + 50) / self.Lf_initial #(self.KP_data.Lf + self.KP_data.CO)
        self.p_1_given_f = self.KP_data.CN / self.Lf_initial
        self.p_0_given_f = 1 - self.p_1_given_f

        # self.p_1_given_s = (np.zeros_like(self.KP_data.DN)) / self.Ls_initial #(self.KP_data.Ls + self.KP_data.DO)
        self.p_1_given_s = self.KP_data.DN / self.Ls_initial
        self.p_0_given_s = 1 - self.p_1_given_s

        # Joint Distributions
        self.p_0_and_f = self.joint_probability(self.p_0_given_f, self.p_f)
        self.p_0_and_s = self.joint_probability(self.p_0_given_s, self.p_s)

        self.p_1_and_f = self.joint_probability(self.p_1_given_f, self.p_f)
        self.p_1_and_s = self.joint_probability(self.p_1_given_s, self.p_s)

        # Marginal Distribution
        self.p_0 = self.marginal_probability(self.p_0_and_f, self.p_0_and_s)
        self.p_1 = self.marginal_probability(self.p_1_and_f, self.p_1_and_s)

    @staticmethod
    def joint_probability(p_o_given_i, p_i):
        p_o_and_i = p_o_given_i * p_i
        return p_o_and_i

    @staticmethod
    def marginal_probability(p_o_and_i1, p_o_and_i2):
        p_o = p_o_and_i1 + p_o_and_i2
        return p_o

    def compute_mutual_information(self):
        term_s_0 = self.p_0_and_s * np.log2(self.p_0_given_s/self.p_0)
        term_s_1 = self.p_1_and_s * np.log2(self.p_1_given_s/self.p_1)
        term_f_0 = self.p_0_and_f * np.log2(self.p_0_given_f/self.p_0)
        term_f_1 = self.p_1_and_f * np.log2(self.p_1_given_f/self.p_1)

        I = term_s_0 + term_s_1 + term_f_0 + term_f_1
        if self.save:
            np.savetxt("mutual_information", I)
        return I

    def plot_I(self):
        plt.plot(self.KP_data.time, self.compute_mutual_information(), label="I(o;i)", linestyle='-', marker='o')
        if self.save:
            plt.savefig("I.pdf", format='pdf')
            plt.close()


class MutualInformationAlternate(MutualInformation):
    def __init__(self, kp_data, save=False):
        MutualInformation.__init__(self, kp_data, save=save)

    def compute_output_entropy(self):
        H_o = -np.nan_to_num((self.p_0 * np.log2(self.p_0) + self.p_1 * np.log2(self.p_1)))
        if self.save:
            np.savetxt("p_0", self.p_0)
            np.savetxt("p_1", self.p_1)
            np.savetxt("output_entropy", H_o)
        return H_o

    def compute_conditional_entropy(self):
        term_s_0 = np.nan_to_num(self.p_0_and_s * np.log2(self.p_0_given_s))
        term_s_1 = np.nan_to_num(self.p_1_and_s * np.log2(self.p_1_given_s))
        term_f_0 = np.nan_to_num(self.p_0_and_f * np.log2(self.p_0_given_f))
        term_f_1 = np.nan_to_num(self.p_1_and_f * np.log2(self.p_1_given_f))

        H_o_given_i = term_s_0 + term_s_1 + term_f_0 + term_f_1
        return H_o_given_i

    def compute_mutual_information(self):
        I = self.compute_output_entropy() + self.compute_conditional_entropy()
        if self.save:
            np.savetxt("mutual_information", I)
        return I

    def plot_I(self, label=None):
        if label:
            plt.plot(self.KP_data.time, self.compute_mutual_information(),
                     label=label, linestyle='-', marker='o')
        else:
            plt.plot(self.KP_data.time, self.compute_mutual_information(),
                     label="{0} steps".format(self.KP_data.num_kp_steps), linestyle='-', marker='o')
        if self.save:
            plt.savefig("I.pdf", format='pdf')
            plt.close()

    def plot_entropies(self):
        plt.plot(self.KP_data.time, self.compute_mutual_information(), label="$I(o;i)$", linestyle='-', marker='o')
        plt.plot(self.KP_data.time, self.compute_output_entropy(), label="$H(o)$", linestyle='-', marker='o')
        plt.plot(self.KP_data.time, -self.compute_conditional_entropy(), label="$H(o|i)$", linestyle='-', marker='o')
        plt.legend()
        if self.save:
            plt.savefig("entropies.pdf", format='pdf')
            plt.close()


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Args for post-processing KP output.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--save', action='store_true', default=True,
                        help='Flag for saving output.')
    parser.add_argument('--steps', dest='steps', action='store', type=int,
                        help="number of KP steps.")
    parser.add_argument('--self', action='store_true', default=False,
                        help='Flag for plotting self')

    args = parser.parse_args()

    # plotkp = PlotKP("mean_traj", args.steps)
    # plotkp.error_fraction()
    # plotkp.save_plot()

    second_order = PlotSecondOrder("mean_traj")
    second_order.save_plot()


    # mutual_information = MutualInformationAlternate(plotkp, save=args.save)
    #
    # print("CN = " + str(mutual_information.KP_data.CN[-5:]))
    # print("DN = " + str(mutual_information.KP_data.DN[-5:]))
    # mutual_information.plot_I()
    # mutual_information.plot_entropies()







