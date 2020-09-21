import argparse
import os

import matplotlib

matplotlib.use('agg')
import matplotlib.pyplot as plt

from compute_ic import InformationCapacity

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot histograms for certain steps",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--steps', dest='steps', action='store', type=int, default=8,
                        help="number of KP steps.")
    parser.add_argument('--fb', dest='fb', action='store', type=float)
    parser.add_argument('--d', dest='d', action='store', type=str)
    parser.add_argument('--xhi', dest='xhi', action='store', type=float)
    parser.add_argument('--xlo', dest='xlo', action='store', type=float)
    parser.add_argument('--lf', dest='lf', action='store', type=int, default=30)

    args = parser.parse_args()
    steps = [int(s) for s in os.path.basename(os.getcwd()) if s.isdigit()][0]
    foreign_directory = [d for d in os.listdir("..") if 'Ls_Lf' in d][0]

    ic_lf = InformationCapacity(foreign_directory=foreign_directory + "/",
                                self_directory="Ls/",
                                limiting="self")
    # print("num_bins = " + str(ic_lf.number_of_bins))
    ic_lf.plot_cn()
    ic_lf.plot_dn()

    plt.xlabel("Concentration of Signaling Output", size=15)
    plt.ylabel("Likelihood", size=20)
    # plt.tick_params(labelsize=15)
    plt.legend(prop={'size': 18})

    if args.xlo and args.xhi:
        plt.xlim(args.xlo, args.xhi)

    steps = 1
    # plt.title("{:.0f} Signaling Steps: C = {:.3f}".format(steps, ic_lf.capacity) + ", $k_{latpp} = 1.0 s^{-1}$",
    #           size=15)
    plt.title("{:.0f} Signaling Step: C = {:.3f}".format(steps, ic_lf.capacity), size=20)
    plt.tight_layout()
    # plt.title("10 Steps: C = {:.3f}".format(ic_lf.calculate_ic()))
    plt.savefig("histograms_{:.0f}.pdf".format(steps), format="pdf")
