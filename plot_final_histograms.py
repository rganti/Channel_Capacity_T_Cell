import argparse

import matplotlib.pyplot as plt

from compute_ic import InformationCapacity

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot histograms for certain steps",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--steps', dest='steps', action='store', type=int, default=8,
                        help="number of KP steps.")
    parser.add_argument('--fb', dest='fb', action='store', type=float)
    parser.add_argument('--d', dest='d', action='store', type=str)

    args = parser.parse_args()
    steps = args.steps
    lf = 30

    # if steps > 8:
    #     if args.fb:
    #         file_path = "{0}_step_k_pos_{1}/".format(steps, args.fb)
    #     else:
    #         file_path = "{0}_step/".format(steps)
    if args.d:
        file_path = args.d
    else:
        file_path = "{0}_step/".format(steps)

    ic_lf = InformationCapacity(foreign_directory=file_path + "Ls_Lf_{0}/".format(lf),
                                self_directory=file_path + "Ls/",
                                limiting="self")
    ic_lf.plot_cn()
    ic_lf.plot_dn()
    plt.legend()

    print("[L_f] = 30")
    plt.xlim(0, 10)

    plt.title("{0} Steps: IC = {1}".format(steps, ic_lf.calculate_ic()))

    plt.savefig(file_path + "histograms.pdf", format="pdf")
