import argparse

from realistic_network import TcrCycleSelfWithForeign, KPRealistic, make_and_cd


class SingleMoleculeKprLs(TcrCycleSelfWithForeign):

    def __init__(self, arguments=None):
        TcrCycleSelfWithForeign.__init__(self, arguments=arguments)

        del self.n_initial['Lf']
        del self.n_initial['Ls']

        # self.p_ligand = [10]

        self.record = []
        self.output = []

    def change_ligand_concentration(self, concentration):
        for i in range(concentration):
            name = "Ls{0}".format(i + 1)
            self.n_initial[name] = 1
            # self.record.append(name)
            self.output.append(name)
            self.k_L_off[name] = self.rate_constants.k_self_off


class SingleMoleculeKprLsLf(SingleMoleculeKprLs):
    def __init__(self, arguments=None):
        SingleMoleculeKprLs.__init__(self, arguments=arguments)

        self.lf = 3

        for i in range(self.lf):
            name = "Lf{0}".format(i + 1)
            self.n_initial[name] = 1
            # self.record.append(name)
            self.output.append(name)
            self.k_L_off[name] = self.rate_constants.k_foreign_off


class KPRealisticSingleMolecule(KPRealistic):
    def __init__(self, self_foreign=False, arguments=None):
        KPRealistic.__init__(self, self_foreign=self_foreign, arguments=arguments)

        if self.self_foreign_flag:
            self.ligand = SingleMoleculeKprLsLf(arguments=arguments)
        else:
            self.ligand = SingleMoleculeKprLs(arguments=arguments)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Submitting job for single step KPR",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--run', action='store_true', default=False,
                        help='Flag for submitting simulations.')
    parser.add_argument('--ss', action='store_true', default=False,
                        help="flag for checking if sims approach steady-state.")
    parser.add_argument('--steps', dest='steps', action='store', type=int, default=0,
                        help="number of KP steps.")
    parser.add_argument('--ls', action='store_true', default=False,
                        help="flag for submitting Ls calculations.")
    parser.add_argument('--ls_lf', dest='ls_lf', action='store', type=int, default=5,
                        help="number of foreign ligands.")
    args = parser.parse_args()

    directory_name = "{0}_step".format(args.steps)
    make_and_cd(directory_name)

    if args.ls:
        sub_directory = "Ls"
        make_and_cd(sub_directory)
        kp = KPRealisticSingleMolecule(arguments=args)

    elif args.ls_lf:
        sub_directory = "Ls_Lf"
        make_and_cd(sub_directory)
        kp = KPRealisticSingleMolecule(self_foreign=True, arguments=args)
    else:
        raise Exception("Need to specify Ls or Ls_Lf")

    kp.main_script(run=args.run)
