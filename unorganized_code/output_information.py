import argparse
import pickle

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from compute_ic import InformationCapacity


class OneDHistogramInformation(object):
    def __init__(self, directory="."):
        self.df = pd.read_csv(directory + "/file_paths", sep='\t', index_col=0)
        self.ic = []
        self.kp_lat_1 = []
        self.file_paths = [directory + "/" + path for path in self.df['file_path'][2:]]
        self.lf = 30

    def compute_variables(self):
        count = 0
        for file_path in self.file_paths:
            ic = InformationCapacity(foreign_directory=file_path + "/Ls_Lf_{0}/".format(self.lf),
                                     self_directory=file_path + "/Ls/", limiting="self")
            C = ic.capacity
            np.savetxt(file_path + "/IC", [C], fmt="%f")
            parameters = pickle.load(open(file_path + "/Ls/parameters.pickle", "rb"))

            self.ic.append(C)
            self.kp_lat_1.append(parameters['kp_on_lat1'])
            count += 1


class HistogramInformation(object):
    def __init__(self, directory=".", steps=8, lf=1):
        self.df = pd.read_csv(directory + "/file_paths", sep='\t', index_col=0)
        self.ic = []
        self.rate_variable = []
        self.on_rate = []

        self.file_paths = [directory + "/" + path for path in self.df['file_path']]
        self.lf = lf
        self.steps = steps

        self.x_array = np.array([1.0, 5.0, 10.0, 50.0, 100.0])
        self.y_array = np.linspace(0.5, 5.0, 10)

    def plot_histograms(self):
        for file_path in self.file_paths:
            print("file_path = " + str(file_path))
            ic = InformationCapacity(foreign_directory=file_path + "/Ls_Lf_{0}/".format(self.lf),
                                     self_directory=file_path + "/Ls/", limiting="self")
            C = ic.capacity
            np.savetxt(file_path + "/IC", [C], fmt="%f")

            ic.plot_cn()
            ic.plot_dn()
            plt.legend()

            print("[L_f] = {0}".format(self.lf))

            plt.title("{:.0f} Steps: C = {:.3f}, P_0_integral = {:.3f}".format(self.steps, C,
                                                                               ic.p0_integral))

            plt.savefig(file_path + "/histograms_lf_{0}.pdf".format(self.lf), format="pdf")
            plt.close()

    def compute_variables(self, rate='kp_on_lat2'):
        k_on = []
        capacity = []
        rate_parameter = []
        count = 0
        for file_path in self.file_paths:
            ic = InformationCapacity(foreign_directory=file_path + "/Ls_Lf_{0}/".format(self.lf),
                                     self_directory=file_path + "/Ls/", limiting="self")
            C = ic.capacity
            np.savetxt(file_path + "/IC", [C], fmt="%f")
            parameters = pickle.load(open(file_path + "/Ls/parameters.pickle", "rb"))

            # ic.plot_cn()
            # ic.plot_dn()
            # plt.legend()
            # plt.xlim(0, 10)
            # plt.title("IC = {0}".format(ic.calculate_ic()))
            #
            # plt.savefig(file_path + "/histograms.pdf", format="pdf")
            # plt.close()

            # lf_mean, lf_iqr = ic.cn_mean_iqr()
            # self.lf_mean.append(lf_mean)
            # self.lf_iqr.append(lf_iqr)
            #
            # ls_mean, ls_iqr = ic.dn_mean_iqr()
            # self.ls_mean.append(ls_mean)
            # self.ls_iqr.append(ls_iqr)

            capacity.append(C)
            rate_parameter.append(parameters[rate])
            k_on.append(parameters['kp_on_lat1'])

            count += 1
            if count == 5:
                self.ic.append(capacity)
                self.rate_variable.append(rate_parameter)
                self.on_rate.append(k_on)

                k_on = []
                rate_parameter = []
                capacity = []
                count = 0

    # def modify_df(self):
    #     self.df['ls_mean'] = self.ls_mean
    #     self.df['lf_mean'] = self.lf_mean
    #
    #     self.df['ls_iqr'] = self.ls_iqr
    #     self.df['lf_iqr'] = self.lf_iqr
    #
    #     self.df['C'] = self.ic
    #
    #     self.df.to_csv("./output_info", sep='\t', float_format='%.3f')

    def main(self):
        self.compute_variables()

        plt.imshow(np.array(self.ic), extent=[1.0, 100.0, 5.0, 0.5], aspect=2)
        plt.xticks(self.x_array)
        plt.yticks(self.y_array)

        plt.ylabel("$k_{on}$")
        plt.xlabel("reduction")

        plt.colorbar(orientation='vertical')

        plt.savefig("capacity_on_rates.pdf", format="pdf")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot histograms for certain steps",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--steps', dest='steps', action='store', type=int, default=8,
                        help="number of KP steps.")
    args = parser.parse_args()

    hist_info = HistogramInformation(steps=args.steps)

    hist_info.plot_histograms()
