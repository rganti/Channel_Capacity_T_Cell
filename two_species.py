#!/usr/bin/python
import argparse
import os
import subprocess
import time

import numpy as np


class SharedCommands(object):

    def __init__(self, n_initial, record, script_name):
        self.n_initial = n_initial
        self.record = record
        self.script_name = script_name

    def define_region(self, f):
        f.write('region World box width 1 height 1 depth 1\n')
        f.write('subvolume edge 1\n\n')

    def initialize(self, f):
        for key, value in self.n_initial.items():
            f.write("new {0} at {1}\n".format(key, value))

    def record_species(self, f):
        for item in self.record:
            f.write("record {0}\n".format(item))
        f.write("\n")


def compile_script(script_name):
    subprocess.call(["ssc", "--save-expanded=network", "{0}".format(script_name)])


class TwoSpecies(object):
    def __init__(self):
        self.k_AB = 0.1
        self.k_BA = 0.2
        self.species = {"A": [self.k_AB, "B"], "B": [self.k_BA, "A"]}
        self.n_initial = {"A": 800, "B": 200}
        self.record = ["A", "B"]
        self.script_name = "two_species_A_B.rxn"
        self.shared = SharedCommands(self.n_initial, self.record, self.script_name)

        self.simulation_name = "two_species_A_B"
        self.num_files = 50
        self.time_step = 0.1

    def define_rxns(self, f):
        for key, value in self.species.items():
            f.write("rxn x:{0} at {1} -> destroy x; new {2}\n".format(key, value[0], value[1]))
        f.write("\n")

    def generate_script(self):
        f = open(self.script_name, "w")

        self.shared.define_region(f)
        self.define_rxns(f)
        self.shared.initialize(f)
        self.shared.record_species(f)

        f.close()

    def generate_qsub(self):
        q = open("qsub.sh", "w")
        q.write("#PBS -m ae\n")
        q.write("#PBS -q short\n")
        q.write("#PBS -V\n")
        q.write("#PBS -l walltime=00:01:00,nodes=1:ppn=1 -N {0}\n\n".format(self.simulation_name))
        q.write("cd $PBS_O_WORKDIR\n\n")
        q.write("EXE_FILE={0}\n".format(self.simulation_name))
        q.write("RUN_TIME={0}\n".format(self.num_files))
        q.write("STEP={0}\n\n".format(self.time_step))
        q.write("for j in {1.." + str(self.num_files) + "}\n")
        q.write("do\n")
        q.write("\t ./$EXE_FILE -e $RUN_TIME -t $STEP > traj_$j\n")
        q.write("done\n")


class TwoWay(object):
    def __init__(self):
        self.dictionary = {}

    def add(self, key, value):
        self.dictionary[key] = value
        self.dictionary[value] = key


class KineticProofreading(object):
    def __init__(self):
        self.num_kp_steps = 1

        self.kappa = 0.0052
        self.k_off_foreign = 0.2
        self.k_off_self = 2.0

        self.n_initial = {"R": 1000, "Lf": 500, "Ls": 500}
        self.record = ["Lf", "Ls", "C0", "D0"]

        self.simulation_name = "kinetic_proofreading"

        self.num_files = 100
        self.run_time = 100
        # self.time_step = 0.1
        self.time_step = 100

        self.forward_rates = {"RLf": self.kappa, "RLs": self.kappa}
        self.reverse_rates = {"C0": self.k_off_foreign, "D0": self.k_off_self}

        self.forward_rxns = [[["R", "Lf"], ["C0"]], [["R", "Ls"], ["D0"]]]
        self.reverse_rxns = [[["C0"], ["R", "Lf"]], [["D0"], ["R", "Ls"]]]

    @staticmethod
    def define_reactions(f, rxn, rate):
        for i in range(len(rxn)):
            input_string = "rxn "
            destroy_string = ""
            input = rxn[i][0]
            output = rxn[i][1]

            rate_key = ""
            for item in input:
                input_string += "{0}:{1} ".format(item.lower(), item)
                destroy_string += "destroy {0}; ".format(item.lower())
                rate_key += item
            output_string = ""
            for item in output:
                output_string += "new {0}; ".format(item)

            rate_string = "at {0} -> ".format(rate[rate_key])

            f.write(input_string + rate_string + destroy_string + output_string + "\n")

    def generate_ssc_script(self):
        script_name = self.simulation_name + ".rxn"
        shared = SharedCommands(self.n_initial, self.record, self.simulation_name + ".rxn")

        f = open(script_name, "w")

        shared.define_region(f)
        self.define_reactions(f, self.forward_rxns, self.forward_rates)
        self.define_reactions(f, self.reverse_rxns, self.reverse_rates)
        f.write("\n")
        shared.initialize(f)
        f.write("\n")
        shared.record_species(f)

        f.close()

    def generate_qsub(self):
        q = open("qsub.sh", "w")
        q.write("#PBS -m ae\n")
        q.write("#PBS -q short\n")
        q.write("#PBS -V\n")
        q.write("#PBS -l walltime=00:01:00,nodes=1:ppn=1 -N {0}\n\n".format(self.simulation_name))
        q.write("cd $PBS_O_WORKDIR\n\n")
        q.write("EXE_FILE={0}\n".format(self.simulation_name))
        q.write("RUN_TIME={0}\n".format(self.run_time))
        q.write("STEP={0}\n\n".format(self.time_step))
        q.write("for j in {1.." + str(self.num_files) + "}\n")
        q.write("do\n")
        # q.write("\t ./$EXE_FILE -e $RUN_TIME -t $STEP > traj_$j\n")
        q.write("\t ./$EXE_FILE -e $RUN_TIME > traj_$j\n")
        q.write("done\n\n")
        q.write("python ~/SSC_python_modules/post_process.py --num_files {0} "
                "--run_time {1} --time_step {2}\n".format(self.num_files, self.run_time, self.time_step))
        # q.write("python ~/SSC_python_modules/plot.py --steps {0}\n".format(self.num_kp_steps))
        q.close()


class ForeignLigand(KineticProofreading):
    def __init__(self):
        KineticProofreading.__init__(self)

        self.inputs = ["R", "Lf"]
        self.n_initial = {"R": 10000, "Lf": 0}
        self.record = ["Lf", "C0"]
        self.simulation_name = "kp_single_species"

        self.forward_rates = {"RLf": self.kappa}
        self.reverse_rates = {"C0": self.k_off_foreign}

        self.forward_rxns = [[["R", "Lf"], ["C0"]]]
        self.reverse_rxns = [[["C0"], ["R", "Lf"]]]

        self.output_symbol = "C"

        self.mu = 20
        self.sigma = 0.5

        self.p_ligand = [int(i) for i in np.round(np.random.normal(self.mu, self.sigma, 1000))]

    def change_ligand_concentration(self, concentration):
        self.n_initial["Lf"] = concentration


class SelfLigand(KineticProofreading):
    def __init__(self):
        KineticProofreading.__init__(self)

        self.inputs = ["R", "Ls"]
        self.n_initial = {"R": 10000, "Ls": 0}
        self.record = ["Ls", "D0"]
        self.simulation_name = "kp_single_species"

        self.forward_rates = {"RLs": self.kappa}
        self.reverse_rates = {"D0": self.k_off_self}

        self.forward_rxns = [[["R", "Ls"], ["D0"]]]
        self.reverse_rxns = [[["D0"], ["R", "Ls"]]]

        self.output_symbol = "D"

        self.mu = 6
        self.sigma = 1.0

        self.p_ligand = [int(i) for i in np.round(np.random.lognormal(self.mu, self.sigma, 1000))]

    def change_ligand_concentration(self, concentration):
        self.n_initial["Ls"] = concentration


class KPSingleSpecies(object):
    def __init__(self, foreign=False):
        if foreign:
            self.ligand = ForeignLigand()
        else:
            self.ligand = SelfLigand()
        self.home_directory = os.getcwd()

        self.k_p = 0.05

        self.num_kp_steps = self.ligand.num_kp_steps
        self.forward_rates = self.ligand.forward_rates
        self.forward_rxns = self.ligand.forward_rxns

        self.reverse_rates = self.ligand.reverse_rates
        self.reverse_rxns = self.ligand.reverse_rxns

        self.symbol = self.ligand.output_symbol
        self.inputs = self.ligand.inputs

        self.record = self.ligand.record

        self.num_samples = 1000

    def generate(self):
        self.ligand.generate_ssc_script()
        compile_script(self.ligand.simulation_name + ".rxn")
        self.ligand.generate_qsub()

    def add_step(self):
        self.num_kp_steps += 1
        self.simulation_name = "kp_steps_" + str(self.num_kp_steps)
        print("Num KP steps = " + str(self.num_kp_steps))

        self.forward_rates[self.symbol + "{0}".format(self.num_kp_steps - 2)] = self.k_p
        self.forward_rxns.append([[self.symbol + "{0}".format(self.num_kp_steps - 2)],
                                  [self.symbol + "{0}".format(self.num_kp_steps - 1)]])

        self.reverse_rates[self.symbol + "{0}".format(self.num_kp_steps - 1)] = self.reverse_rates[self.symbol + "0"]
        self.reverse_rxns.append([[self.symbol + "{0}".format(self.num_kp_steps - 1)], self.inputs])

        self.record.append(self.symbol + "{0}".format(self.num_kp_steps - 1))

    def main_script(self, run=False):
        sample = []
        for i in range(self.num_samples):
            directory = "sample_" + str(i)
            s = self.ligand.p_ligand[i]
            sample.append(s)
            self.ligand.change_ligand_concentration(s)
            self.ligand.simulation_name = "kp_single_species_" + str(i)

            os.makedirs(directory)
            print("Made " + directory)
            os.chdir(directory)
            print("Changed into directory: " + str(os.getcwd()))

            self.generate()
            if run:
                (stdout, stderr) = subprocess.Popen(["qsub {0}".format("qsub.sh")], shell=True, stdout=subprocess.PIPE,
                                                    cwd=os.getcwd()).communicate()
            os.chdir(self.home_directory)
            time.sleep(0.1)

        np.savetxt("Ligand_concentrations", sample, fmt='%f')
        np.savetxt("Ligand_concentrations_sorted", np.sort(sample), fmt='%f')


# class KPAddSteps(object):
#     def __init__(self, foreign=False):
#         self.KPScheme = KPSingleSpecies(foreign=foreign)
#         self.k_p = 0.05
#
#         self.num_kp_steps = self.KPScheme.ligand.num_kp_steps
#         self.forward_rates = self.KPScheme.ligand.forward_rates
#         self.forward_rxns = self.KPScheme.ligand.forward_rxns
#
#         self.reverse_rates = self.KPScheme.ligand.reverse_rates
#         self.reverse_rxns = self.KPScheme.ligand.reverse_rxns
#
#         self.symbol = self.KPScheme.ligand.output_symbol
#         self.inputs = self.KPScheme.ligand.inputs
#
#         self.record = self.KPScheme.ligand.record
#
#     def add_step(self):
#         self.num_kp_steps += 1
#         self.simulation_name = "kp_steps_" + str(self.num_kp_steps)
#         print("Num KP steps = " + str(self.num_kp_steps))
#
#         self.forward_rates[self.symbol + "{0}".format(self.num_kp_steps - 2)] = self.k_p
#         self.forward_rxns.append([[self.symbol + "{0}".format(self.num_kp_steps - 2)],
#                                   [self.symbol + "{0}".format(self.num_kp_steps - 1)]])
#
#         self.reverse_rates[self.symbol + "{0}".format(self.num_kp_steps - 1)] = self.reverse_rates[self.symbol + "0"]
#         self.reverse_rxns.append([[self.symbol + "{0}".format(self.num_kp_steps - 1)], self.inputs])
#
#         self.record.append(self.symbol + "{0}".format(self.num_kp_steps - 1))


class KPMultiStep(KineticProofreading):
    def __init__(self):
        KineticProofreading.__init__(self)
        self.k_p = 0.05

    def add_step(self):
        self.num_kp_steps += 1
        self.simulation_name = "kp_steps_" + str(self.num_kp_steps)
        print("Num KP steps = " + str(self.num_kp_steps))
        for i in ["C", "D"]:
            self.forward_rates[i + "{0}".format(self.num_kp_steps - 2)] = self.k_p
            self.forward_rxns.append([[i + "{0}".format(self.num_kp_steps - 2)], [i + "{0}".format(self.num_kp_steps - 1)]])
            if i == "C":
                self.reverse_rates[i + "{0}".format(self.num_kp_steps - 1)] = self.k_off_foreign
                self.reverse_rxns.append([[i + "{0}".format(self.num_kp_steps - 1)], ["R", "Lf"]])
            elif i == "D":
                self.reverse_rates[i + "{0}".format(self.num_kp_steps - 1)] = self.k_off_self
                self.reverse_rxns.append([[i + "{0}".format(self.num_kp_steps - 1)], ["R", "Ls"]])
            self.record.append(i + "{0}".format(self.num_kp_steps - 1))

    def negative_feedback_loop(self):
        if self.num_kp_steps == 3:
            k_negative = 1.0
            self.reverse_rxns += [[['C1', 'C0'], ['R', 'Lf']]]
            self.reverse_rates['C1C0'] = k_negative

            self.reverse_rxns += [[['D1', 'D0'], ['R', 'Ls']]]
            self.reverse_rates['D1D0'] = k_negative

    def positive_feedback_loop(self):
        if self.num_kp_steps == 3:
            self.forward_rxns += [[['C1', 'C2'], ['C2']]]
            self.forward_rates['C1C2'] = 0.01

            self.forward_rxns += [[['D1', 'D2'], ['D2']]]
            self.forward_rates['D1D2'] = 0.01


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Submitting job for calculating P(C0) as function of steps",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--run', action='store_true', default=False,
                        help='Flag for submitting simulations.')
    args = parser.parse_args()

    if "foreign" in os.getcwd():
        kp = KPSingleSpecies(foreign=True)
    else:
        kp = KPSingleSpecies()
    kp.add_step()
    kp.add_step()
    kp.main_script(run=args.run)


    # kp_multistep = KPMultiStep()
    # # kp_multistep.add_step()
    # # kp_multistep.add_step()
    # kp_multistep.generate_ssc_script()
    # print("Number of steps = " + str(kp_multistep.num_kp_steps))
    #
    # compile_script(kp_multistep.simulation_name + ".rxn")
    # kp_multistep.generate_qsub()
