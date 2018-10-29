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
    parser.add_argument('--xhi', dest='xhi', action='store', type=float)
    parser.add_argument('--xlo', dest='xlo', action='store', type=float)
    parser.add_argument('--lf', dest='lf', action='store', type=int, default=30)

    args = parser.parse_args()
    steps = args.steps
    lf = args.lf

    # if steps > 8:
    #     if args.fb:
    #         file_path = "{0}_step_k_pos_{1}/".format(steps, args.fb)
    #     else:
    #         file_path = "{0}_step/".format(steps)
    # if args.d:
    #     file_path = args.d
    # else:
    #     file_path = "{0}_step/".format(steps)

    ic_lf = InformationCapacity(foreign_directory="Ls_Lf_{0}/".format(lf),
                                self_directory="Ls/",
                                limiting="self")
    ic_lf.plot_cn()
    ic_lf.plot_dn()
    plt.legend()

    print("[L_f] = {0}".format(lf))

    if args.xlo and args.xhi:
        plt.xlim(args.xlo, args.xhi)

    plt.title("{:.0f} Steps: C = {:.3f}, P_0_integral = {:.3f}".format(steps, ic_lf.calculate_ic(), ic_lf.p_0_integral))

    plt.savefig("histograms_lf_{0}.pdf".format(lf), format="pdf")
