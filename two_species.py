#!/usr/bin/python
import argparse
import datetime
import os
import subprocess

import numpy as np

from simulation_parameters import DefineRegion


class SharedCommands(object):

    def __init__(self, n_initial, record):
        self.n_initial = n_initial
        self.record = record

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
        self.shared = SharedCommands(self.n_initial, self.record)

        self.regions = DefineRegion()

        self.simulation_name = "two_species_A_B"
        self.num_files = 50
        self.time_step = 0.1

    def define_rxns(self, f):
        for key, value in self.species.items():
            f.write("rxn x:{0} at {1} -> destroy x; new {2}\n".format(key, value[0], value[1]))
        f.write("\n")

    def generate_script(self):
        f = open(self.script_name, "w")

        self.regions.define_region(f)
        self.define_rxns(f)
        self.shared.initialize(f)
        self.shared.record_species(f)

        f.close()

    def generate_qsub(self):
        q = open("qsub.sh", "w")
        q.write("#PBS -m ae\n")
        q.write("#PBS -q short\n")
        q.write("#PBS -V\n")
        q.write("#PBS -l walltime=00:02:00,nodes=1:ppn=1 -N {0}\n\n".format(self.simulation_name))
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


class KPBindingParameters(object):
    def __init__(self):
        self.k_tcr_on = 0.0052
        self.k_foreign_off = 0.2
        self.k_self_off = 10.0 * self.k_foreign_off

        self.k_p = 0.05
        self.k_p_off_foreign = 0.5
        self.k_p_off_self = 10.0 * self.k_p_off_foreign


class SelfWithForeign(object):
    def __init__(self, arguments=None):
        self.arguments = arguments
        self.rate_constants = KPBindingParameters()

        self.n_initial = {"R": 10000, "Lf": 20, "Ls": 0}
        self.record = ["Lf", "Ls", "C0", "D0"]

        self.simulation_name = "kp_competing"

        self.forward_rates = {"RLf": self.rate_constants.k_tcr_on, "RLs": self.rate_constants.k_tcr_on}
        self.reverse_rates = {"C0": self.rate_constants.k_foreign_off, "D0": self.rate_constants.k_self_off}

        self.forward_rxns = [[["R", "Lf"], ["C0"]], [["R", "Ls"], ["D0"]]]
        self.reverse_rxns = [[["C0"], ["R", "Lf"]], [["D0"], ["R", "Ls"]]]

        self.num_kp_steps = 1
        self.num_samples = 0
        if self.arguments:
            if self.arguments.test or self.arguments.ss:
                self.num_samples = 1
            else:
                self.num_samples = 1000

        self.mu = 6
        self.sigma = 1.0
        self.p_ligand = [int(i) for i in np.round(np.random.lognormal(self.mu, self.sigma, self.num_samples))]

        self.output = ["C", "D"]

    def change_ligand_concentration(self, concentration):
        self.n_initial["Ls"] = concentration

    def modify_forward_reverse(self, reactants, products, forward_rate, reverse_rate):
        self.forward_rates[''.join(reactants)] = forward_rate
        self.forward_rxns.append([reactants, products])

        self.reverse_rates[''.join(products)] = reverse_rate
        self.reverse_rxns.append([products, reactants])

        self.record.append(''.join(products))

    def increment_step(self):
        self.num_kp_steps += 1
        self.simulation_name = "kp_steps_" + str(self.num_kp_steps)
        print("Num KP steps = " + str(self.num_kp_steps))

    def add_step_1(self):
        self.increment_step()
        for i in self.output:
            if i == "C":
                k_p_off = self.rate_constants.k_p_off_foreign
            elif i == "D":
                k_p_off = self.rate_constants.k_p_off_self

            self.modify_forward_reverse([i + "0"], [i + "1"], self.rate_constants.k_p, k_p_off)

    def add_step_2(self):
        self.increment_step()
        for i in self.output:
            if i == "C":
                k_p_off = self.rate_constants.k_p_off_foreign
            elif i == "D":
                k_p_off = self.rate_constants.k_p_off_self

            self.modify_forward_reverse([i + "1"], [i + "2"], self.rate_constants.k_p, k_p_off)

    def add_step_3(self):
        self.increment_step()
        for i in self.output:
            if i == "C":
                k_p_off = self.rate_constants.k_p_off_foreign
            elif i == "D":
                k_p_off = self.rate_constants.k_p_off_self

            self.modify_forward_reverse([i + "2"], [i + "3"], self.rate_constants.k_p, k_p_off)


class ReversibleSelfLigand(SelfWithForeign):
    def __init__(self):
        SelfWithForeign.__init__(self)
        del self.n_initial['Lf']

        self.n_initial = {"R": 10000, "Ls": 0}
        self.record = ["Ls", "D0"]
        self.simulation_name = "kp_ls"

        self.forward_rates = {"RLs": self.rate_constants.k_tcr_on}
        self.reverse_rates = {"D0": self.rate_constants.k_self_off}

        self.forward_rxns = [[["R", "Ls"], ["D0"]]]
        self.reverse_rxns = [[["D0"], ["R", "Ls"]]]

        self.output = ["D"]


class ForeignLigand(object):
    def __init__(self, arguments=None):
        self.arguments = arguments
        self.rate_constants = KPBindingParameters()

        self.inputs = ["R", "Lf"]
        self.n_initial = {"R": 10000, "Lf": 0}
        self.record = ["Lf", "C0"]
        self.simulation_name = "kp_lf"

        self.forward_rates = {"RLf": self.rate_constants.k_tcr_on}
        self.reverse_rates = {"C0": self.rate_constants.k_foreign_off}

        self.forward_rxns = [[["R", "Lf"], ["C0"]]]
        self.reverse_rxns = [[["C0"], ["R", "Lf"]]]

        self.symbol = "C"

        self.num_kp_steps = 1
        self.num_samples = 0
        if self.arguments:
            if self.arguments.test or self.arguments.ss:
                self.num_samples = 5
            else:
                self.num_samples = 1000

        self.mu = 20
        self.sigma = 0.5

        self.p_ligand = [int(i) for i in np.round(np.random.normal(self.mu, self.sigma, self.num_samples))]

    def change_ligand_concentration(self, concentration):
        self.n_initial["Lf"] = concentration


class SelfLigand(object):
    def __init__(self, arguments=None):
        self.arguments = arguments
        self.rate_constants = KPBindingParameters()

        self.inputs = ["R", "Ls"]
        self.n_initial = {"R": 10000, "Ls": 0}
        self.record = ["Ls", "D0"]
        self.simulation_name = "kp_ls"

        self.forward_rates = {"RLs": self.rate_constants.k_tcr_on}
        self.reverse_rates = {"D0": self.rate_constants.k_self_off}

        self.forward_rxns = [[["R", "Ls"], ["D0"]]]
        self.reverse_rxns = [[["D0"], ["R", "Ls"]]]

        self.symbol = "D"

        self.num_kp_steps = 1
        if self.arguments:
            if self.arguments.test or self.arguments.ss:
                self.num_samples = 5
            else:
                self.num_samples = 1000

        self.mu = 6
        self.sigma = 1.0

        self.p_ligand = [int(i) for i in np.round(np.random.lognormal(self.mu, self.sigma, self.num_samples))]

    def change_ligand_concentration(self, concentration):
        self.n_initial["Ls"] = concentration


def add_to_network(n, reactants, products, rate):
    n.write(" + ".join(reactants))
    n.write(" -> {0} ".format(rate))
    n.write(" + ".join(products))
    n.write("\n")


class KPSingleSpecies(object):
    def __init__(self, foreign=False, self_foreign=False, arguments=None):
        self.foreign_flag = foreign
        self.self_foreign_flag = self_foreign
        self.arguments = arguments

        self.regions = DefineRegion()

        if self.foreign_flag:
            self.ligand = ForeignLigand()
        elif self.self_foreign_flag:
            self.ligand = SelfWithForeign()
        else:
            self.ligand = ReversibleSelfLigand()

        self.num_files = 100
        self.run_time = 100
        self.simulation_time = 2
        self.single_molecule = False

        self.home_directory = os.getcwd()

        self.num_kp_steps = 1

    def set_simulation_time(self):
        simulation_time = self.run_time * (20.0 / 1000)

        return simulation_time

    def set_time_step(self):
        if self.arguments:
            if self.arguments.ss:
                time_step = 1.0
            else:
                time_step = self.run_time
        else:
            time_step = self.run_time
        return time_step

    @staticmethod
    def define_reactions(f, rxn, rate, n):
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

            rate_key += "_"
            output_string = ""
            for item in output:
                output_string += "new {0}; ".format(item)
                rate_key += item

            rate_string = "at {0} -> ".format(rate[rate_key])

            f.write(input_string + rate_string + destroy_string + output_string + "\n")
            add_to_network(n, input, output, rate[rate_key])

    def generate_ssc_script(self, simulation_name):
        script_name = simulation_name + ".rxn"
        shared = SharedCommands(self.ligand.n_initial, self.ligand.record)

        f = open(script_name, "w")
        n = open("ordered_network", "w")

        self.regions.define_region(f)
        f.write("-- Forward reactions \n")
        n.write("# Forward Reactions \n")
        self.define_reactions(f, self.ligand.forward_rxns, self.ligand.forward_rates, n)

        n.write("\n# Reverse Reactions \n")
        f.write("\n-- Reverse reactions \n")
        self.define_reactions(f, self.ligand.reverse_rxns, self.ligand.reverse_rates, n)
        f.write("\n")
        shared.initialize(f)
        f.write("\n")
        shared.record_species(f)

        n.close()
        f.close()

    def generate_qsub(self, simulation_name, time_step):
        q = open("qsub.sh", "w")
        q.write("#PBS -m ae\n")
        q.write("#PBS -q short\n")
        q.write("#PBS -V\n")
        q.write("#PBS -l walltime={1},nodes=1:ppn=1 -N {0}\n\n".format(simulation_name,
                                                                       datetime.timedelta(
                                                                           minutes=self.set_simulation_time())))
        q.write("cd $PBS_O_WORKDIR\n\n")
        q.write("echo $PBS_JOBID > job_id\n")
        q.write("EXE_FILE={0}\n".format(simulation_name))
        q.write("RUN_TIME={0}\n".format(self.run_time))
        q.write("STEP={0}\n\n".format(time_step))
        q.write("for j in {1.." + str(self.num_files) + "}\n")
        q.write("do\n")
        if time_step == self.run_time:
            q.write("\t ./$EXE_FILE -e $RUN_TIME > traj_$j\n")
        else:
            q.write("\t ./$EXE_FILE -e $RUN_TIME -t $STEP > traj_$j\n")
        q.write("done\n\n")
        q.write("python ~/SSC_python_modules/post_process.py --num_files {0} "
                "--run_time $RUN_TIME --time_step $STEP\n".format(self.num_files))
        if self.single_molecule:
            q.write("wait \n")
            q.write("python ~/SSC_python_modules/kp_sm_post_process.py \n")
        if self.arguments.ss:
            q.write("python ~/SSC_python_modules/plot.py \n")
        q.close()

    def generate(self, simulation_name, time_step):
        self.generate_ssc_script(simulation_name)
        compile_script(simulation_name + ".rxn")
        self.generate_qsub(simulation_name, time_step)

    def single_add_step(self):
        self.num_kp_steps += 1
        self.ligand.simulation_name = "kp_steps_" + str(self.num_kp_steps)
        print("Num KP steps = " + str(self.num_kp_steps))

        self.ligand.forward_rates[
            self.ligand.symbol + "{0}".format(self.num_kp_steps - 2)] = self.ligand.rate_constants.k_p
        self.ligand.forward_rxns.append([[self.ligand.symbol + "{0}".format(self.num_kp_steps - 2)],
                                  [self.ligand.symbol + "{0}".format(self.num_kp_steps - 1)]])

        self.ligand.reverse_rates[self.ligand.symbol + "{0}".format(self.num_kp_steps - 1)] = self.ligand.reverse_rates[
            self.ligand.symbol + "0"]
        self.ligand.reverse_rxns.append(
            [[self.ligand.symbol + "{0}".format(self.num_kp_steps - 1)], self.ligand.inputs])

        self.ligand.record.append(self.ligand.symbol + "{0}".format(self.num_kp_steps - 1))

    def competing_add_step(self):
        self.num_kp_steps += 1
        self.ligand.simulation_name = "kp_steps_" + str(self.num_kp_steps)
        print("Num KP steps = " + str(self.num_kp_steps))
        for i in ["C", "D"]:
            self.ligand.forward_rates[i + "{0}".format(self.num_kp_steps - 2)] = self.ligand.rate_constants.k_p
            self.ligand.forward_rxns.append([[i + "{0}".format(self.num_kp_steps - 2)],
                                      [i + "{0}".format(self.num_kp_steps - 1)]])
            if i == "C":
                self.ligand.reverse_rates[
                    i + "{0}".format(self.num_kp_steps - 1)] = self.ligand.rate_constants.k_foreign_off
                self.ligand.reverse_rxns.append([[i + "{0}".format(self.num_kp_steps - 1)], ["R", "Lf"]])
            elif i == "D":
                self.ligand.reverse_rates[
                    i + "{0}".format(self.num_kp_steps - 1)] = self.ligand.rate_constants.k_self_off
                self.ligand.reverse_rxns.append([[i + "{0}".format(self.num_kp_steps - 1)], ["R", "Ls"]])
            self.ligand.record.append(i + "{0}".format(self.num_kp_steps - 1))

    def add_step(self):
        if self.self_foreign_flag:
            self.competing_add_step()
        else:
            self.single_add_step()

    def main_script(self, run=False):
        sample = []
        for i in range(self.ligand.num_samples):
            directory = "sample_" + str(i)
            s = self.ligand.p_ligand[i]
            sample.append(s)
            self.ligand.change_ligand_concentration(s)
            simulation_name = self.ligand.simulation_name + "_" + str(i)

            os.makedirs(directory)
            print("Made " + directory)
            os.chdir(directory)
            print("Changed into directory: " + str(os.getcwd()))

            if self.ligand.num_kp_steps > 6:
                self.run_time = 1000

            self.generate(simulation_name, self.set_time_step())
            if run:
                (stdout, stderr) = subprocess.Popen(["qsub {0}".format("qsub.sh")], shell=True, stdout=subprocess.PIPE,
                                                    cwd=os.getcwd()).communicate()
            os.chdir(self.home_directory)

        np.savetxt("Ligand_concentrations", sample, fmt='%f')
        np.savetxt("Ligand_concentrations_sorted", np.sort(sample), fmt='%f')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Submitting job for calculating P(C0) as function of steps",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--run', action='store_true', default=False,
                        help='Flag for submitting simulations.')
    parser.add_argument('--test', action='store_true', default=False,
                        help="flag for testing.")
    parser.add_argument('--ss', action='store_true', default=False,
                        help="flag for checking if sims approach steady-state.")
    args = parser.parse_args()

    if "Ls_Lf" in os.getcwd():
        kp = KPSingleSpecies(self_foreign=True, arguments=args)
    elif "Lf" in os.getcwd():
        kp = KPSingleSpecies(foreign=True, arguments=args)
    elif "Ls" in os.getcwd():
        kp = KPSingleSpecies()
    else:
        raise Exception("Incorrect Directory labeling. Specify (Ls, Lf, Ls_Lf)")

    # kp.add_step()
    # kp.add_step()

    kp.ligand.add_step_1()
    # kp.ligand.add_step_2()
    # kp.ligand.add_step_3()
    # kp.ligand.add_step_4()
    # kp.ligand.add_step_5()
    # kp.ligand.add_step_6()
    # kp.ligand.add_step_7()
    # kp.ligand.add_step_8()
    # kp.ligand.add_step_9()

    kp.main_script(run=args.run)
