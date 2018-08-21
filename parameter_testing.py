import argparse
import os
import subprocess
import time

import numpy as np
import pandas as pd

from post_process import load
from realistic_network import TcrCycleSelfLigand, make_and_cd, TcrCycleSelfWithForeign
from two_species import KPSingleSpecies


def lognuniform(low=0, high=1, size=None, base=10):
    return np.power(base, np.random.uniform(low, high, size))


class KPParameterTest(KPSingleSpecies):
    def __init__(self, self_foreign=False, arguments=None):
        KPSingleSpecies.__init__(self, self_foreign=self_foreign, arguments=arguments)

        if self.self_foreign_flag:
            print("initializing Self and Foreign")
            self.ligand = TcrCycleSelfWithForeign(arguments=arguments)
        else:
            print("initializing Self Only")
            self.ligand = TcrCycleSelfLigand(arguments=arguments)

        self.ligand.n_initial["Ls"] = 800
        if self.arguments.run:
            self.run_time = self.arguments.run
        else:
            self.run_time = 100

        self.forward_rates = self.ligand.forward_rates
        self.forward_rxns = self.ligand.forward_rxns

        self.reverse_rates = self.ligand.reverse_rates
        self.reverse_rxns = self.ligand.reverse_rxns

        self.record = self.ligand.record

        self.simulation_time = self.set_simulation_time()

    def set_simulation_time(self):
        simulation_time = self.run_time * (12.0 / 1000)

        return simulation_time


class ParameterTesting(object):
    def __init__(self, arguments=None):
        self.arguments = arguments
        self.parameter_list = np.sort(lognuniform(low=-4, high=2, size=20))

        if self.arguments.ss or self.arguments.test:
            x = self.parameter_list[:1]
            y = self.parameter_list[-1:]
            self.parameter_list = [x[0], y[0]]

        self.ligand_directories = ["Ls", "Ls_Lf_{0}".format(self.arguments.ls_lf)]

        self.file_list = []

        self.foreign_file_list = []
        self.self_file_list = []

    def wait_for_simulations(self):
        while True:
            list1 = []

            for file_path in self.file_list:
                list1.append(os.path.isfile(file_path + "mean_traj"))
                print(str(list1))

            if all(list1):
                # All elements are True. Therefore all the files exist. Run %run commands
                break
            else:
                # At least one element is False. Therefore not all the files exist. Run FTP commands again
                time.sleep(5)  # wait 30 seconds before checking again

    def simple_kp(self, run=False):
        sub_directory = os.getcwd()

        for ligand_directory in self.ligand_directories:
            if ligand_directory == "Ls":
                kp_ligand = KPParameterTest(arguments=self.arguments)
                self.self_file_list.append("{0}/".format(ligand_directory))
            else:
                kp_ligand = KPParameterTest(self_foreign=True, arguments=self.arguments)
                self.foreign_file_list.append("{0}/".format(ligand_directory))

            self.file_list = self.foreign_file_list + self.self_file_list
            make_and_cd(ligand_directory)

            kp_ligand.ligand.add_step_0()
            self.add_steps(kp_ligand)

            kp_ligand.generate(kp_ligand.ligand.simulation_name, kp_ligand.set_time_step())

            if run:
                (stdout, stderr) = subprocess.Popen(["qsub {0}".format("qsub.sh")], shell=True, stdout=subprocess.PIPE,
                                                    cwd=os.getcwd()).communicate()

            os.chdir(sub_directory)

        if run:
            self.wait_for_simulations()
            self.compute_hopfield_error()

    @staticmethod
    def change_parameter(kp_ligand, parameter):
        kp_ligand.ligand.rate_constants.k_lck_off_zap_R = parameter

    def parameter_testing(self, run=False):
        home_directory = os.getcwd()
        for i in range(len(self.parameter_list)):
            dir_name = "parameter_{0}".format(i)
            make_and_cd(dir_name)
            sub_directory = os.getcwd()

            for ligand_directory in self.ligand_directories:
                if ligand_directory == "Ls":
                    kp_ligand = KPParameterTest(arguments=self.arguments)
                    self.self_file_list.append("{0}/{1}/".format(dir_name, ligand_directory))
                else:
                    kp_ligand = KPParameterTest(self_foreign=True, arguments=self.arguments)
                    self.foreign_file_list.append("{0}/{1}/".format(dir_name, ligand_directory))

                self.file_list = self.foreign_file_list + self.self_file_list
                make_and_cd(ligand_directory)

                self.change_parameter(kp_ligand, self.parameter_list[i])
                # kp_ligand.ligand.rate_constants.k_lck_off_zap_R = self.parameter_list[i]
                kp_ligand.ligand.add_step_0()

                self.add_steps(kp_ligand)

                kp_ligand.generate(kp_ligand.ligand.simulation_name, kp_ligand.set_time_step())

                if run:
                    (stdout, stderr) = subprocess.Popen(["qsub {0}".format("qsub.sh")], shell=True,
                                                        stdout=subprocess.PIPE,
                                                        cwd=os.getcwd()).communicate()

                os.chdir(sub_directory)

            os.chdir(home_directory)

        np.savetxt("rates", self.parameter_list, fmt='%f')

        if run:
            self.wait_for_simulations()
            self.compute_hopfield_error(parameter_test=True)

    def check_columns(self, file_path):
        column_names = load(file_path + "column_names")[0].split()
        if "Ls_Lf_{0}".format(self.arguments.ls_lf) in self.ligand_directories:
            if "Lf" in column_names[-1] or "Ls" in column_names[-1]:
                self_foreign = True
            else:
                self_foreign = False
        else:
            self_foreign = False
        return self_foreign

    def compute_hopfield_error(self, parameter_test=False):
        eta_array = []
        self_output = []
        foreign_output = []

        ligand_output = self.check_columns(self.foreign_file_list[0])
        for i in range(len(self.self_file_list)):
            self_trajectory = np.loadtxt(self.self_file_list[i] + "mean_traj")
            self_output.append(self_trajectory[-1])

            foreign_trajectory = np.loadtxt(self.foreign_file_list[i] + "mean_traj")
            if ligand_output:
                eta_array.append(self_trajectory[-1] / (foreign_trajectory[-1] + foreign_trajectory[-2]))
                foreign_output.append(foreign_trajectory[-1] + foreign_trajectory[-2])
            else:
                eta_array.append(self_trajectory[-1] / foreign_trajectory[-1])
                foreign_output.append(foreign_trajectory[-1])

        print("All arrays processed.")
        if parameter_test:
            d = {"rates": self.parameter_list, "eta": eta_array, "Lf_output": foreign_output, "Ls_output": self_output}
        else:
            d = {"eta": eta_array, "Lf_output": foreign_output, "Ls_output": self_output}
        df = pd.DataFrame(data=d)

        df.to_csv("eta", "\t", float_format='%.6f')

    def process_output(self):
        # for i in range(len(self.parameter_list)):
        #     dir_name = "parameter_{0}".format(i)
        for ligand_directory in self.ligand_directories:
            if ligand_directory == "Ls":
                self.self_file_list.append("{0}/".format(ligand_directory))
            else:
                self.foreign_file_list.append("{0}/".format(ligand_directory))

        self.compute_hopfield_error()

    def add_steps(self, kp_ligand):
        if self.arguments.steps > 0:
            kp_ligand.ligand.add_cycle(kp_ligand.ligand.cycle_1)
        if self.arguments.steps > 1:
            kp_ligand.ligand.add_cycle(kp_ligand.ligand.cycle_2)
        if self.arguments.steps > 2:
            kp_ligand.ligand.add_cycle(kp_ligand.ligand.cycle_3)
        if self.arguments.steps > 3:
            kp_ligand.ligand.add_cycle(kp_ligand.ligand.cycle_4)
        if self.arguments.steps > 4:
            kp_ligand.ligand.add_negative_feedback()
        if self.arguments.steps > 5:
            kp_ligand.ligand.add_cycle(kp_ligand.ligand.cycle_6)
        if self.arguments.steps > 6:
            kp_ligand.ligand.add_step_7()
        if self.arguments.steps > 7:
            kp_ligand.ligand.add_step_8()
        if self.arguments.steps > 8:
            kp_ligand.ligand.add_positive_feedback()


class ParameterTestingSecondOrder(ParameterTesting):
    def __init__(self, arguments=None):
        ParameterTesting.__init__(self, arguments=arguments)
        self.parameter_list = [float(i) for i in np.sort(lognuniform(low=-4, high=-1, size=10))]

        if self.arguments.ss or self.arguments.test:
            x = self.parameter_list[:1]
            y = self.parameter_list[-1:]
            self.parameter_list = [x[0], y[0]]

    @staticmethod
    def change_parameter(kp_ligand, parameter):
        kp_ligand.ligand.rate_constants.k_product_on = parameter


class ParameterTestingSecondOrderNumbers(ParameterTestingSecondOrder):
    def __init__(self, arguments=None):
        ParameterTestingSecondOrder.__init__(self, arguments=arguments)
        self.parameter_list = [int(i) for i in np.sort(lognuniform(low=1, high=4, size=10))]

        if self.arguments.ss or self.arguments.test:
            x = self.parameter_list[:1]
            y = self.parameter_list[-1:]
            self.parameter_list = [x[0], y[0]]

    @staticmethod
    def change_parameter(kp_ligand, parameter):
        kp_ligand.ligand.n_initial["LAT"] = parameter


def make_and_cd_p_test(directory_name):
    count = 1
    while os.path.exists(directory_name):
        new_directory_name = "old_" + directory_name + "_" + str(count)
        if os.path.exists(new_directory_name):
            count += 1
            continue
        else:
            os.rename(directory_name, new_directory_name)
            print("Changed " + directory_name + " to " + new_directory_name)
            break

    os.mkdir(directory_name)
    os.chdir(directory_name)
    print("Changed into directory: " + str(os.getcwd()))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Submitting job for calculating P(C0) as function of steps",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--run', dest='run', action='store', type=int,
                        help='Flag for submitting simulations.')
    parser.add_argument('--p_test', action='store_true', default=False,
                        help='Flag for submitting parameter testing simulations.')
    parser.add_argument('--ss', dest='ss', action='store_true', default=False,
                        help="flag for checking if sims approach steady-state.")
    parser.add_argument('--test', action='store_true', default=False,
                        help="flag for testing.")
    parser.add_argument('--steps', dest='steps', action='store', type=int, default=3,
                        help="number of KP steps.")
    parser.add_argument('--ls_lf', dest='ls_lf', action='store', type=int, default=30,
                        help="number of foreign ligands.")

    args = parser.parse_args()
    # print(str(args.steps))

    kp_parameter_testing = ParameterTestingSecondOrder(arguments=args)
    # kp_parameter_testing.process_output()

    if args.p_test:
        directory_name = "{0}_step_p_test".format(args.steps)
        make_and_cd_p_test(directory_name)

        kp_parameter_testing.parameter_testing(run=args.run)
    else:
        if args.ss:
            directory_name = "{0}_step_ss".format(args.steps)
        else:
            directory_name = "{0}_step".format(args.steps)
        make_and_cd_p_test(directory_name)

        kp_parameter_testing.simple_kp(run=args.run)
