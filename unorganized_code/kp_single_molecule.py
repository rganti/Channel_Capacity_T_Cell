import argparse

from realistic_network import TcrCycleSelfWithForeign, KPRealistic, make_and_cd
from simulation_parameters import MembraneBindingParameters
from ssc_tcr_membrane import MembraneSharedCommands


class SingleMoleculeKprLs(TcrCycleSelfWithForeign):

    def __init__(self, arguments=None):
        TcrCycleSelfWithForeign.__init__(self, arguments=arguments)
        self.rate_constants = MembraneBindingParameters()

        del self.n_initial['Lf']
        del self.n_initial['Ls']

        self.record = []
        self.output = []

        self.diffusion_flag = True

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

        self.single_molecule = True

    def define_diffusion(self, f):

        for key in self.ligand.diffusion_rate_dict.keys():
            if self.ligand.diffusion_loc_dict[key] == "Plasma":
                f.write("diffusion {0} at {1} in {2}\n".format(key, self.ligand.diffusion_rate_dict[key],
                                                               self.ligand.diffusion_loc_dict[key]))
            else:
                f.write(
                    "diffusion {0} at {1} in {2}, Cytosol<->Plasma\n".format(key, self.ligand.diffusion_rate_dict[key],
                                                                             self.ligand.diffusion_loc_dict[key]))

    def generate_ssc_script(self, simulation_name):
        script_name = simulation_name + ".rxn"
        shared = MembraneSharedCommands(self.ligand.n_initial, self.ligand.record, self.ligand.diffusion_loc_dict)

        f = open(script_name, "w")
        n = open("ordered_network", "w")

        self.regions.define_membrane_region(f)
        f.write("-- Forward reactions \n")
        n.write("# Forward Reactions \n")
        self.define_reactions(f, self.ligand.forward_rxns, self.ligand.forward_rates, n)

        n.write("\n# Reverse Reactions \n")
        f.write("\n-- Reverse reactions \n")
        self.define_reactions(f, self.ligand.reverse_rxns, self.ligand.reverse_rates, n)
        f.write("\n")

        if self.ligand.diffusion_flag:
            f.write("\n-- Diffusion \n")
            self.define_diffusion(f)
            f.write("\n")

        shared.initialize(f)
        f.write("\n")
        shared.record_species(f)

        n.close()
        f.close()


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
    parser.add_argument('--ls_lf', dest='ls_lf', action='store', type=int, default=3,
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
        # make_and_cd(sub_directory)
        kp = KPRealisticSingleMolecule(self_foreign=True, arguments=args)
    else:
        raise Exception("Need to specify Ls or Ls_Lf")

    kp.main_script(run=args.run)
