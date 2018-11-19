import argparse
import os

from realistic_network import KPRealistic, TcrCycleSelfWithForeign, make_and_cd
from simulation_parameters import MembraneInitialConcentrations, MembraneBindingParameters
from two_species import SharedCommands


class MembraneSharedCommands(SharedCommands):
    def __init__(self, n_initial, record, species_location):
        SharedCommands.__init__(self, n_initial, record)
        self.species_location = species_location

    def initialize(self, f):
        for key, value in self.n_initial.items():
            f.write("new {0} at {1} in {2}\n".format(key, value, self.species_location[key]))


class MembraneTcrSelfWithForeign(TcrCycleSelfWithForeign):
    def __init__(self, arguments=None):
        TcrCycleSelfWithForeign.__init__(self, arguments=arguments)

        self.initial = MembraneInitialConcentrations()
        self.rate_constants = MembraneBindingParameters()

        self.n_initial = {"R": self.initial.r_0, "Lf": self.arguments.ls_lf, "Ls": 1000}

        if self.arguments.steps > 0:
            self.n_initial["Lck"] = self.initial.lck_0

        if self.arguments.steps > 2:
            self.n_initial["Zap"] = self.initial.zap_0

        if self.arguments.steps > 4:
            self.n_initial["LAT"] = self.initial.lat_0

        self.diffusion_flag = True
        self.num_samples = 3
        self.p_ligand = [50, 100, 150]


class MembraneTcrSelfLigand(MembraneTcrSelfWithForeign):
    def __init__(self, arguments=None):
        MembraneTcrSelfWithForeign.__init__(self, arguments=arguments)

        del self.n_initial['Lf']
        self.n_initial["Ls"] = 1000
        self.record = ["Ls"]
        self.output = ["Ls"]


class KpMembraneSpecies(KPRealistic):
    def __init__(self, self_foreign=False, arguments=None):
        KPRealistic.__init__(self, self_foreign=self_foreign, arguments=arguments)

        if self.self_foreign_flag:
            self.ligand = MembraneTcrSelfWithForeign(arguments=arguments)
        else:
            self.ligand = MembraneTcrSelfLigand(arguments=arguments)

        self.run_time = 500

        self.forward_rates = self.ligand.forward_rates
        self.forward_rxns = self.ligand.forward_rxns

        self.reverse_rates = self.ligand.reverse_rates
        self.reverse_rxns = self.ligand.reverse_rxns

        self.record = self.ligand.record

    def define_diffusion(self, f):

        for key in self.ligand.diffusion_rate_dict.keys():
            if self.ligand.diffusion_loc_dict[key] == "Plasma":
                f.write("diffusion {0} at {1} in {2}\n".format(key, self.ligand.diffusion_rate_dict[key],
                                                               self.ligand.diffusion_loc_dict[key]))

    def generate_ssc_script(self, simulation_name):
        script_name = simulation_name + ".rxn"
        shared = MembraneSharedCommands(self.ligand.n_initial, self.record, self.ligand.diffusion_loc_dict)

        f = open(script_name, "w")
        n = open("ordered_network", "w")

        self.regions.define_membrane_region(f)
        f.write("-- Forward reactions \n")
        n.write("# Forward Reactions \n")
        self.define_reactions(f, self.forward_rxns, self.forward_rates, n)

        n.write("\n# Reverse Reactions \n")
        f.write("\n-- Reverse reactions \n")
        self.define_reactions(f, self.reverse_rxns, self.reverse_rates, n)
        f.write("\n")

        f.write("\n-- Diffusion \n")
        self.define_diffusion(f)
        f.write("\n")

        shared.initialize(f)
        f.write("\n")
        shared.record_species(f)

        n.close()
        f.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Submitting job for calculating P(C0) as function of steps",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--run', action='store_true', default=False,
                        help='Flag for submitting simulations.')
    parser.add_argument('--ss', action='store_true', default=False,
                        help="flag for checking if sims approach steady-state.")
    parser.add_argument('--steps', dest='steps', action='store', type=int, default=0,
                        help="number of KP steps.")
    parser.add_argument('--ls_lf', dest='ls_lf', action='store', type=int, default=5,
                        help="number of foreign ligands.")
    args = parser.parse_args()

    directory_name = "{0}_step_spatial".format(args.steps)
    make_and_cd(directory_name)

    sub_directories = ["Ls"]  # , "Ls_Lf_{0}".format(args.ls_lf)]
    home_directory = os.getcwd()

    for sub_directory in sub_directories:
        # make_and_cd(sub_directory)
        if sub_directory == "Ls":
            kp = KpMembraneSpecies(arguments=args)
        else:
            kp = KpMembraneSpecies(self_foreign=True, arguments=args)

        kp.ligand.add_step_0()

        if args.steps > 0:
            kp.ligand.add_cycle(kp.ligand.cycle_1)
        if args.steps > 1:
            kp.ligand.add_cycle(kp.ligand.cycle_2)

        # if args.steps > 2:
        #     kp.ligand.add_cycle(kp.ligand.cycle_3)
        # if args.steps > 3:
        #     kp.ligand.add_cycle(kp.ligand.cycle_4)
        # if args.steps > 5:
        #     kp.ligand.add_cycle(kp.ligand.cycle_6)
        # if args.steps > 6:
        #     kp.ligand.add_cycle(kp.ligand.cycle_7)
        # if args.steps > 7:
        #     kp.ligand.add_cycle(kp.ligand.cycle_8)

        kp.main_script(run=args.run)
        os.chdir(home_directory)
