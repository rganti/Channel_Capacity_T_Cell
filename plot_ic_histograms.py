import argparse

import matplotlib.pyplot as plt
import numpy as np

from compute_ic import InformationCapacity

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot histograms for certain steps",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--steps', dest='steps', action='store', type=int, default=8,
                        help="number of KP steps.")
    parser.add_argument('--fb', dest='fb', action='store', type=float)

    args = parser.parse_args()
    steps = args.steps
    lf = 30

    num_steps = []
    capacity = []
    for i in range(1, 10):
        if i == 5:
            continue

        if i > 8:
            if args.fb:
                file_path = "{0}_step_k_pos_{1}/".format(i, args.fb)
            else:
                file_path = "{0}_step/".format(i)
        else:
            file_path = "{0}_step/".format(i)

        num_steps.append(i)
        ic_lf = InformationCapacity(foreign_directory=file_path + "Ls_Lf_{0}/".format(lf),
                                    self_directory=file_path + "Ls/", limiting="self")
        ic_lf.plot_cn()
        ic_lf.plot_dn()

        if i == 1:
            xhi = 1500
        elif i == 2:
            xhi = 1000
        elif i == 3:
            xhi = 600
        elif i == 4:
            xhi = 400
        elif i > 8:
            xhi = 30
        else:
            xhi = 10

        plt.xlim(0, xhi)
        plt.title("{0} Steps: IC = {1}".format(i, ic_lf.calculate_ic()))
        plt.legend()
        plt.savefig(file_path + "{0}_step.pdf".format(i), format="pdf")
        plt.close()
        #
        capacity.append(ic_lf.calculate_ic())

        plt.plot(num_steps, capacity, linestyle='-', marker='o', label="[Lf] = 30")

        plt.legend()
        plt.xlabel("Number of Steps", size=15)
        plt.ylabel("C (bits)", size=15)
        plt.xlim(0, 10)
        plt.ylim(0, 1)
        plt.locator_params(axis='x', nbins=11)
        plt.savefig(file_path + "ic_{0}_step.pdf".format(i), format="pdf")
        plt.close()

    np.savetxt("num_steps", num_steps, fmt='%f')
    np.savetxt("ic", capacity, fmt='%f')
