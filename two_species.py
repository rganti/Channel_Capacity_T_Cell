#!/usr/bin/python
import argparse
import os
import subprocess

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


class BindingParameters(object):
    def __init__(self):
        self.k_tcr_on = 0.0022
        self.k_foreign_off = 0.01
        self.k_self_off = 10.0 * self.k_foreign_off

        self.k_lck_on = 0.297
        self.k_lck_off = 10.0

        self.k_p_on = 5.0
        self.k_p_off = 1.0

        self.k_zap_on = 135.0
        self.k_zap_off = 0.11

        self.k_lat_on = 0.1
        self.k_lat_off = 0.1

        self.k_plc_on = 161.2
        self.k_plc_off = 1.0

        self.k_p_plc_on = 0.01

        self.k_pip_on = 0.1
        self.k_pip_off = 0.2


class TcrSelfWithForeign(object):
    def __init__(self):
        self.rate_constants = BindingParameters()

        self.n_initial = {"R": 10000, "Lck": 10000, "Zap": 2300000, "Lat": 3000, "Plc": 460000, "Pip": 45000,
                          "Lf": 20, "Ls": 0}
        self.record = ["Lf", "Ls", "C_RL", "D_RL"]

        self.simulation_name = "kp_competing"

        self.forward_rates = {"RLf": self.rate_constants.k_tcr_on, "RLs": self.rate_constants.k_tcr_on}
        self.reverse_rates = {"C_RL": self.rate_constants.k_foreign_off, "D_RL": self.rate_constants.k_self_off}

        self.forward_rxns = [[["R", "Lf"], ["C_RL"]], [["R", "Ls"], ["D_RL"]]]
        self.reverse_rxns = [[["C_RL"], ["R", "Lf"]], [["D_RL"], ["R", "Ls"]]]

        self.mu = 6
        self.sigma = 1.0
        if args.test or args.ss:
            self.num_samples = 5
        else:
            self.num_samples = 1000

        self.p_ligand = [int(i) for i in np.round(np.random.lognormal(self.mu, self.sigma, self.num_samples))]

        self.num_kp_steps = 1
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
            self.modify_forward_reverse([i + "_RL", "Lck"], [i + "_RL1"], self.rate_constants.k_lck_on,
                                        self.rate_constants.k_lck_off)

    def add_step_2(self):
        self.increment_step()
        for i in self.output:
            self.modify_forward_reverse([i + "_RL1"], [i + "_RL2"], self.rate_constants.k_p_on,
                                        self.rate_constants.k_p_off)

    def add_step_3(self):
        self.increment_step()
        for i in self.output:
            self.modify_forward_reverse([i + "_RL2", "Zap"], [i + "_Zap"], self.rate_constants.k_zap_on,
                                        self.rate_constants.k_zap_off)

    def add_step_4(self):
        self.increment_step()
        for i in self.output:
            self.modify_forward_reverse([i + "_Zap"], [i + "_Zap1"], self.rate_constants.k_p_on,
                                        self.rate_constants.k_p_off)

    def add_step_5(self):
        self.increment_step()
        for i in self.output:
            self.modify_forward_reverse([i + "_Zap1", "Lat"], [i + "_Lat"], self.rate_constants.k_lat_on,
                                        self.rate_constants.k_lat_off)

    def add_step_6(self):
        self.increment_step()
        for i in self.output:
            self.modify_forward_reverse([i + "_Lat"], [i + "_Lat1"], self.rate_constants.k_p_on,
                                        self.rate_constants.k_p_off)

    def add_step_7(self):
        self.increment_step()
        for i in self.output:
            self.modify_forward_reverse([i + "_Lat1", "Plc"], [i + "_Plc"], self.rate_constants.k_plc_on,
                                        self.rate_constants.k_plc_off)

    def add_step_8(self):
        self.increment_step()
        for i in self.output:
            self.modify_forward_reverse([i + "_Plc"], [i + "_Plc1"], self.rate_constants.k_p_plc_on,
                                        self.rate_constants.k_p_off)

    def add_step_9(self):
        self.increment_step()
        for i in self.output:
            self.modify_forward_reverse([i + "_Plc1", "Pip"], [i + "_Pip"], self.rate_constants.k_pip_on,
                                        self.rate_constants.k_pip_off)

    # def add_step_1(self):
    #     self.num_kp_steps += 1
    #     self.simulation_name = "kp_steps_" + str(self.num_kp_steps)
    #     print("Num KP steps = " + str(self.num_kp_steps))
    #     for i in self.output:
    #         self.forward_rates[i + "_RLLck"] = self.rate_constants.k_lck_on
    #         self.forward_rxns.append([[i + "_RL", "Lck"], [i + "_RL1"]])
    #
    #         self.reverse_rates[i + "_RL1"] = self.rate_constants.k_lck_off
    #         self.reverse_rxns.append([[i + "_RL1"], [i + "_RL", "Lck"]])
    #
    #         self.record.append(i + "_RL1")

    # def add_step_2(self):
    #     self.num_kp_steps += 1
    #     self.simulation_name = "kp_steps_" + str(self.num_kp_steps)
    #     print("Num KP steps = " + str(self.num_kp_steps))
    #     for i in self.output:
    #         self.forward_rates[i + "_RL1"] = self.rate_constants.k_p_on
    #         self.forward_rxns.append([[i + "_RL1"], [i + "_RL2"]])
    #
    #         self.reverse_rates[i + "_RL2"] = self.rate_constants.k_p_off
    #         self.reverse_rxns.append([[i + "_RL2"], [i + "_RL1"]])
    #
    #         self.record.append(i + "_RL2")

    # def add_step_3(self):
    #     self.num_kp_steps += 1
    #     self.simulation_name = "kp_steps_" + str(self.num_kp_steps)
    #     print("Num KP steps = " + str(self.num_kp_steps))
    #     for i in self.output:
    #         self.forward_rates[i + "_RL2Zap"] = self.rate_constants.k_zap_on
    #         self.forward_rxns.append([[i + "_RL2", "Zap"], [i + "_Zap"]])
    #
    #         self.reverse_rates[i + "_Zap"] = self.rate_constants.k_zap_off
    #         self.reverse_rxns.append([[i + "_Zap"], [i + "_RL2", "Zap"]])
    #
    #         self.record.append(i + "_Zap")

    # def add_step_4(self):
    #     self.num_kp_steps += 1
    #     self.simulation_name = "kp_steps_" + str(self.num_kp_steps)
    #     print("Num KP steps = " + str(self.num_kp_steps))
    #     for i in self.output:
    #         self.forward_rates[i + "_Zap"] = self.rate_constants.k_p_on
    #         self.forward_rxns.append([[i + "_Zap"], [i + "_Zap1"]])
    #
    #         self.reverse_rates[i + "_Zap1"] = self.rate_constants.k_p_off
    #         self.reverse_rxns.append([[i + "_Zap1"], [i + "_Zap"]])
    #
    #         self.record.append(i + "_Zap1")

    # def add_step_5(self):
    #     self.num_kp_steps += 1
    #     self.simulation_name = "kp_steps_" + str(self.num_kp_steps)
    #     print("Num KP steps = " + str(self.num_kp_steps))
    #     for i in self.output:
    #         self.forward_rates[i + "_Zap1"] = self.rate_constants.k_lat_on
    #         self.forward_rxns.append([[i + "_Zap1", "Lat"], [i + "_Lat"]])
    #
    #         self.reverse_rates[i + "_Lat"] = self.rate_constants.k_lat_off
    #         self.reverse_rxns.append([[i + "_Lat"], [i + "_Zap1", "Lat"]])
    #
    #         self.record.append(i + "_Lat")


class TcrSelfLigand(TcrSelfWithForeign):
    def __init__(self):
        TcrSelfWithForeign.__init__(self)

        del self.n_initial['Lf']
        self.record = ["Ls", "D_RL"]
        self.simulation_name = "kp_ls"

        self.forward_rates = {"RLs": self.rate_constants.k_tcr_on}
        self.reverse_rates = {"D_RL": self.rate_constants.k_self_off}

        self.forward_rxns = [[["R", "Ls"], ["D_RL"]]]
        self.reverse_rxns = [[["D_RL"], ["R", "Ls"]]]

        self.output = ["D"]


class KPBindingParameters(object):
    def __init__(self):
        self.k_tcr_on = 0.0052
        self.k_foreign_off = 0.2
        # self.k_self_off = 10.0 * self.k_foreign_off
        self.k_self_off = self.k_foreign_off
        self.k_p = 0.05


class SelfWithForeign(object):
    def __init__(self):
        self.rate_constants = KPBindingParameters()

        self.n_initial = {"R": 10000, "Lf": 20, "Ls": 0}
        self.record = ["Lf", "Ls", "C0", "D0"]

        self.simulation_name = "kp_competing"

        self.forward_rates = {"RLf": self.rate_constants.k_tcr_on, "RLs": self.rate_constants.k_tcr_on}
        self.reverse_rates = {"C0": self.rate_constants.k_foreign_off, "D0": self.rate_constants.k_self_off}

        self.forward_rxns = [[["R", "Lf"], ["C0"]], [["R", "Ls"], ["D0"]]]
        self.reverse_rxns = [[["C0"], ["R", "Lf"]], [["D0"], ["R", "Ls"]]]

        self.mu = 6
        self.sigma = 1.0

        self.num_kp_steps = 1
        if args.test or args.ss:
            self.num_samples = 5
        else:
            self.num_samples = 1000

        self.p_ligand = [int(i) for i in np.round(np.random.lognormal(self.mu, self.sigma, 1000))]

    def change_ligand_concentration(self, concentration):
        self.n_initial["Ls"] = concentration


class ForeignLigand(object):
    def __init__(self):
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
        if args.test or args.ss:
            self.num_samples = 5
        else:
            self.num_samples = 1000

        self.mu = 20
        self.sigma = 0.5

        self.p_ligand = [int(i) for i in np.round(np.random.normal(self.mu, self.sigma, 1000))]

    def change_ligand_concentration(self, concentration):
        self.n_initial["Lf"] = concentration


class SelfLigand(object):
    def __init__(self):
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
        if args.test or args.ss:
            self.num_samples = 5
        else:
            self.num_samples = 1000

        self.mu = 6
        self.sigma = 1.0

        self.p_ligand = [int(i) for i in np.round(np.random.lognormal(self.mu, self.sigma, 1000))]

    def change_ligand_concentration(self, concentration):
        self.n_initial["Ls"] = concentration


class KPSingleSpecies(object):
    def __init__(self, foreign=False, self_foreign=False):
        self.foreign_flag = foreign
        self.self_foreign_flag = self_foreign

        if self.foreign_flag:
            self.ligand = ForeignLigand()
        elif self.self_foreign_flag:
            # self.ligand = TcrSelfWithForeign()
            self.ligand = SelfWithForeign()
        else:
            # self.ligand = TcrSelfLigand()
            self.ligand = SelfLigand()

        self.num_files = 100
        self.run_time = 100
        if args.ss:
            self.time_step = 0.1
        else:
            self.time_step = self.run_time

        self.home_directory = os.getcwd()

        self.num_kp_steps = 1
        self.forward_rates = self.ligand.forward_rates
        self.forward_rxns = self.ligand.forward_rxns

        self.reverse_rates = self.ligand.reverse_rates
        self.reverse_rxns = self.ligand.reverse_rxns

        self.record = self.ligand.record

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

    def generate_ssc_script(self, simulation_name):
        script_name = simulation_name + ".rxn"
        shared = SharedCommands(self.ligand.n_initial, self.record, script_name)

        f = open(script_name, "w")

        shared.define_region(f)
        f.write("-- Forward reactions \n")
        self.define_reactions(f, self.forward_rxns, self.forward_rates)
        f.write("\n-- Reverse reactions \n")
        self.define_reactions(f, self.reverse_rxns, self.reverse_rates)
        f.write("\n")
        shared.initialize(f)
        f.write("\n")
        shared.record_species(f)

        f.close()

    def generate_qsub(self, simulation_name):
        q = open("qsub.sh", "w")
        q.write("#PBS -m ae\n")
        q.write("#PBS -q short\n")
        q.write("#PBS -V\n")
        q.write("#PBS -l walltime=00:02:00,nodes=1:ppn=1 -N {0}\n\n".format(simulation_name))
        q.write("cd $PBS_O_WORKDIR\n\n")
        q.write("echo $PBS_JOBID > job_id\n")
        q.write("EXE_FILE={0}\n".format(simulation_name))
        q.write("RUN_TIME={0}\n".format(self.run_time))
        q.write("STEP={0}\n\n".format(self.time_step))
        q.write("for j in {1.." + str(self.num_files) + "}\n")
        q.write("do\n")
        if self.time_step == self.run_time:
            q.write("\t ./$EXE_FILE -e $RUN_TIME > traj_$j\n")
        else:
            q.write("\t ./$EXE_FILE -e $RUN_TIME -t $STEP > traj_$j\n")
        q.write("done\n\n")
        q.write("python ~/SSC_python_modules/post_process.py --num_files {0} "
                "--run_time {1} --time_step {2}\n".format(self.num_files, self.run_time, self.time_step))
        # q.write("python ~/SSC_python_modules/plot.py --steps {0}\n".format(self.num_kp_steps))
        q.close()

    def generate(self, simulation_name):
        self.generate_ssc_script(simulation_name)
        compile_script(simulation_name + ".rxn")
        self.generate_qsub(simulation_name)

    def single_add_step(self):
        self.num_kp_steps += 1
        self.ligand.simulation_name = "kp_steps_" + str(self.num_kp_steps)
        print("Num KP steps = " + str(self.num_kp_steps))

        self.forward_rates[self.ligand.symbol + "{0}".format(self.num_kp_steps - 2)] = self.ligand.rate_constants.k_p
        self.forward_rxns.append([[self.ligand.symbol + "{0}".format(self.num_kp_steps - 2)],
                                  [self.ligand.symbol + "{0}".format(self.num_kp_steps - 1)]])

        self.reverse_rates[self.ligand.symbol + "{0}".format(self.num_kp_steps - 1)] = self.reverse_rates[
            self.ligand.symbol + "0"]
        self.reverse_rxns.append([[self.ligand.symbol + "{0}".format(self.num_kp_steps - 1)], self.ligand.inputs])

        self.record.append(self.ligand.symbol + "{0}".format(self.num_kp_steps - 1))

    def competing_add_step(self):
        self.num_kp_steps += 1
        self.ligand.simulation_name = "kp_steps_" + str(self.num_kp_steps)
        print("Num KP steps = " + str(self.num_kp_steps))
        for i in ["C", "D"]:
            self.forward_rates[i + "{0}".format(self.num_kp_steps - 2)] = self.ligand.rate_constants.k_p
            self.forward_rxns.append([[i + "{0}".format(self.num_kp_steps - 2)],
                                      [i + "{0}".format(self.num_kp_steps - 1)]])
            if i == "C":
                self.reverse_rates[i + "{0}".format(self.num_kp_steps - 1)] = self.ligand.rate_constants.k_foreign_off
                self.reverse_rxns.append([[i + "{0}".format(self.num_kp_steps - 1)], ["R", "Lf"]])
            elif i == "D":
                self.reverse_rates[i + "{0}".format(self.num_kp_steps - 1)] = self.ligand.rate_constants.k_self_off
                self.reverse_rxns.append([[i + "{0}".format(self.num_kp_steps - 1)], ["R", "Ls"]])
            self.record.append(i + "{0}".format(self.num_kp_steps - 1))

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
                if args.ss:
                    self.time_step = 0.1
                else:
                    self.time_step = self.run_time

            self.generate(simulation_name)
            if run:
                (stdout, stderr) = subprocess.Popen(["qsub {0}".format("qsub.sh")], shell=True, stdout=subprocess.PIPE,
                                                    cwd=os.getcwd()).communicate()
            os.chdir(self.home_directory)

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


# class KPMultiStep(SelfWithForeign):
#     def __init__(self):
#         SelfWithForeign.__init__(self)
#         self.k_p = 0.05
#
#     def add_step(self):
#         self.num_kp_steps += 1
#         self.simulation_name = "kp_steps_" + str(self.num_kp_steps)
#         print("Num KP steps = " + str(self.num_kp_steps))
#         for i in ["C", "D"]:
#             self.forward_rates[i + "{0}".format(self.num_kp_steps - 2)] = self.k_p
#             self.forward_rxns.append([[i + "{0}".format(self.num_kp_steps - 2)], [i + "{0}".format(self.num_kp_steps - 1)]])
#             if i == "C":
#                 self.reverse_rates[i + "{0}".format(self.num_kp_steps - 1)] = self.k_off_foreign
#                 self.reverse_rxns.append([[i + "{0}".format(self.num_kp_steps - 1)], ["R", "Lf"]])
#             elif i == "D":
#                 self.reverse_rates[i + "{0}".format(self.num_kp_steps - 1)] = self.k_off_self
#                 self.reverse_rxns.append([[i + "{0}".format(self.num_kp_steps - 1)], ["R", "Ls"]])
#             self.record.append(i + "{0}".format(self.num_kp_steps - 1))
#
#     def negative_feedback_loop(self):
#         if self.num_kp_steps == 3:
#             k_negative = 1.0
#             self.reverse_rxns += [[['C1', 'C0'], ['R', 'Lf']]]
#             self.reverse_rates['C1C0'] = k_negative
#
#             self.reverse_rxns += [[['D1', 'D0'], ['R', 'Ls']]]
#             self.reverse_rates['D1D0'] = k_negative
#
#     def positive_feedback_loop(self):
#         if self.num_kp_steps == 3:
#             self.forward_rxns += [[['C1', 'C2'], ['C2']]]
#             self.forward_rates['C1C2'] = 0.01
#
#             self.forward_rxns += [[['D1', 'D2'], ['D2']]]
#             self.forward_rates['D1D2'] = 0.01


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Submitting job for calculating P(C0) as function of steps",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--run', action='store_true', default=False,
                        help='Flag for submitting simulations.')
    parser.add_argument('--test', action='store_true', default=False,
                        help="flag for testing.")
    parser.add_argument('--ss', action='store_true', default=False,
                        help="falge for checking if sims approach steady-state.")
    args = parser.parse_args()

    if "Ls_Lf" in os.getcwd():
        kp = KPSingleSpecies(self_foreign=True)
    elif "Lf" in os.getcwd():
        kp = KPSingleSpecies(foreign=True)
    elif "Ls" in os.getcwd():
        kp = KPSingleSpecies()
    else:
        raise Exception("Incorrect Directory labeling. Specify (Ls, Lf, Ls_Lf)")

    kp.add_step()
    kp.add_step()

    # kp.ligand.add_step_1()
    # kp.ligand.add_step_2()
    # kp.ligand.add_step_3()
    # kp.ligand.add_step_4()
    # kp.ligand.add_step_5()
    # kp.ligand.add_step_6()
    # kp.ligand.add_step_7()
    # kp.ligand.add_step_8()
    # kp.ligand.add_step_9()

    kp.main_script(run=args.run)
