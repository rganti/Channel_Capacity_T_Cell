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
            C = ic.calculate_ic()
            np.savetxt(file_path + "/IC", [C], fmt="%f")
            parameters = pickle.load(open(file_path + "/Ls/parameters.pickle", "rb"))

            self.ic.append(C)
            self.kp_lat_1.append(parameters['kp_on_lat1'])
            count += 1


class HistogramInformation(object):
    def __init__(self):
        self.df = pd.read_csv("./file_paths", sep='\t', index_col=0)
        self.ic = []
        self.kp_lat_2 = []

        self.file_paths = self.df['file_path'][10:]

        self.lf = 30

    def compute_variables(self):
        k_on = []
        kp_2 = []
        count = 0
        for file_path in self.file_paths:
            ic = InformationCapacity(foreign_directory=file_path + "/Ls_Lf_{0}/".format(self.lf),
                                     self_directory=file_path + "/Ls/", limiting="self")
            C = ic.calculate_ic()
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

            k_on.append(C)
            kp_2.append(parameters['kp_on_lat2'])
            count += 1
            if count == 5:
                self.ic.append(k_on)
                self.kp_lat_2.append(kp_2)
                kp_2 = []
                k_on = []
                count = 0

    def modify_df(self):
        self.df['ls_mean'] = self.ls_mean
        self.df['lf_mean'] = self.lf_mean

        self.df['ls_iqr'] = self.ls_iqr
        self.df['lf_iqr'] = self.lf_iqr

        self.df['C'] = self.ic

        self.df.to_csv("./output_info", sep='\t', float_format='%.3f')

    def main(self):
        self.compute_variables()

        plt.imshow(np.array(self.ic), extent=[1.0, 100.0, 5.0, 0.5], aspect=2)
        plt.xticks(np.array([1.0, 5.0, 10.0, 50.0, 100.0]))
        plt.yticks(np.linspace(0.5, 5.0, 10))

        plt.ylabel("$k_{on}$")
        plt.xlabel("reduction")

        plt.colorbar(orientation='vertical')

        plt.savefig("capacity_on_rates.pdf", format="pdf")

        # self.modify_df()


if __name__ == "__main__":
    hist_info = HistogramInformation()

    hist_info.main()
