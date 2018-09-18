import pandas as pd

from compute_ic import InformationCapacity


class HistogramInformation(object):
    def __init__(self):
        self.df = pd.read_csv("./file_paths", sep='\t', index_col=0)
        self.ic = []
        self.lf_mean = []
        self.ls_mean = []
        self.lf_iqr = []
        self.ls_iqr = []

        self.lf = 30

    def compute_variables(self):
        for file_path in self.df['file_path']:
            ic = InformationCapacity(foreign_directory=file_path + "/Ls_Lf_{0}/".format(self.lf),
                                     self_directory=file_path + "/Ls/", limiting="self")

            # ic.plot_cn()
            # ic.plot_dn()
            # plt.legend()
            # plt.xlim(0, 10)
            # plt.title("IC = {0}".format(ic.calculate_ic()))
            #
            # plt.savefig(file_path + "/histograms.pdf", format="pdf")
            # plt.close()

            lf_mean, lf_iqr = ic.cn_mean_iqr()
            self.lf_mean.append(lf_mean)
            self.lf_iqr.append(lf_iqr)

            ls_mean, ls_iqr = ic.dn_mean_iqr()
            self.ls_mean.append(ls_mean)
            self.ls_iqr.append(ls_iqr)

            self.ic.append(ic.calculate_ic())

    def modify_df(self):
        self.df['ls_mean'] = self.ls_mean
        self.df['lf_mean'] = self.lf_mean

        self.df['ls_iqr'] = self.ls_iqr
        self.df['lf_iqr'] = self.lf_iqr

        self.df['C'] = self.ic

        self.df.to_csv("./output_info", sep='\t', float_format='%.3f')

    def main(self):
        self.compute_variables()
        self.modify_df()


if __name__ == "__main__":
    hist_info = HistogramInformation()

    hist_info.main()
