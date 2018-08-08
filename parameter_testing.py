import argparse
import os
import subprocess
import time

import numpy as np
from matplotlib import pyplot as plt

from post_process import load
from realistic_network import TcrCycleForeignLigand, TcrCycleSelfLigand
from two_species import KPSingleSpecies


def lognuniform(low=0, high=1, size=None, base=10):
    return np.power(base, np.random.uniform(low, high, size))

def make_and_cd(directory_name):
    os.makedirs(directory_name)
    print("Made " + directory_name)
    os.chdir(directory_name)
    print("Changed into directory: " + str(os.getcwd()))


def plot_eta_previous(file_path, eta_class):
    eta = np.loadtxt(file_path + "eta")

    foreign_column = load(file_path + "Lf/column_names")
    foreign_column_names = foreign_column[0].split()

    self_column = load(file_path + "Ls/column_names")
    self_column_names = self_column[0].split()

    plt.plot(eta_class.rates, np.zeros_like(eta_class.rates) + eta,
             label="{0}/{1}".format(self_column_names[-1], foreign_column_names[-1]),
             linestyle='-', marker='o')


class CheckEquality(object):
    def __init__(self, file_path_1, file_path_2):
        self.file_1 = load(file_path_1)
        self.file_2 = load(file_path_2)

    def check(self):
        i = 0
        while i < len(self.file_1):
            x = self.file_1[i] == self.file_2[i]

            if not x:
                print("File 1: " + self.file_1[i])
                print("File 2: " + self.file_2[i])
                break

            i += 1

        return x


class PlotEta(object):
    def __init__(self, file_path="./"):
        self.eta = np.loadtxt(file_path + "eta")
        self.rates = np.loadtxt(file_path + "rates")

        self.foreign_column = load(file_path + "parameter_0/Lf/column_names")
        self.foreign_column_names = self.foreign_column[0].split()

        self.self_column = load(file_path + "parameter_0/Ls/column_names")
        self.self_column_names = self.self_column[0].split()

    def plot_eta(self):
        plt.plot(self.rates, self.eta,
                 label="{0}/{1}".format(self.self_column_names[-1], self.foreign_column_names[-1]),
                 linestyle='-', marker='o')


class KPParameterTest(KPSingleSpecies):
    def __init__(self, self_foreign=False, arguments=None):
        KPSingleSpecies.__init__(self, self_foreign=self_foreign, arguments=arguments)
        # self.run_time = 500

        if self.self_foreign_flag:
            self.ligand = TcrCycleForeignLigand(arguments=arguments)
        else:
            self.ligand = TcrCycleSelfLigand(arguments=arguments)

        self.forward_rates = self.ligand.forward_rates
        self.forward_rxns = self.ligand.forward_rxns

        self.reverse_rates = self.ligand.reverse_rates
        self.reverse_rxns = self.ligand.reverse_rxns

        self.record = self.ligand.record

        self.simulation_time = 5


class ParameterTesting(object):
    def __init__(self, arguments=None):
        self.arguments = arguments
        self.parameter_list = np.sort(lognuniform(low=-4, high=2, size=20))

        if self.arguments.ss or self.arguments.test:
            x = self.parameter_list[:1]
            y = self.parameter_list[-1:]
            self.parameter_list = [x[0], y[0]]

        self.ligand_directories = ["Lf", "Ls"]

        self.file_list = []

        self.foreign_file_list = []
        self.self_file_list = []

    def add_steps(self, kp_ligand):
        if self.arguments.steps > 0:
            kp_ligand.ligand.add_cycle_1_first_order()
        if self.arguments.steps > 1:
            kp_ligand.ligand.add_cycle_2_first_order()
        if self.arguments.steps > 2:
            kp_ligand.ligand.add_cycle_3_first_order()
        if self.arguments.steps > 3:
            kp_ligand.ligand.add_negative_feedback()

    def wait_for_simulations(self):
        while True:
            list1 = []

            for file_path in self.file_list:
                list1.append(os.path.isfile(file_path))
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
            if ligand_directory == "Lf":
                kp_ligand = KPParameterTest(self_foreign=True, arguments=self.arguments)
                self.foreign_file_list.append("{0}/mean_traj".format(ligand_directory))
            else:
                kp_ligand = KPParameterTest(arguments=self.arguments)
                self.self_file_list.append("{0}/mean_traj".format(ligand_directory))

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
            compute_output(self.foreign_file_list)
            compute_output(self.self_file_list)

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
                if ligand_directory == "Lf":
                    kp_ligand = KPParameterTest(self_foreign=True, arguments=self.arguments)
                    self.foreign_file_list.append("{0}/{1}/mean_traj".format(dir_name, ligand_directory))
                else:
                    kp_ligand = KPParameterTest(arguments=self.arguments)
                    self.self_file_list.append("{0}/{1}/mean_traj".format(dir_name, ligand_directory))

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
            self.compute_hopfield_error()
            compute_output(self.foreign_file_list)
            compute_output(self.self_file_list)

    def compute_hopfield_error(self):
        eta_array = []
        for i in range(len(self.self_file_list)):
            self_trajectory = np.loadtxt(self.self_file_list[i])
            foreign_trajectory = np.loadtxt(self.foreign_file_list[i])

            eta_array.append(self_trajectory[-1] / foreign_trajectory[-1])

        np.savetxt("eta", eta_array, fmt='%f')

    def process_output(self):
        for i in range(len(self.parameter_list)):
            dir_name = "parameter_{0}".format(i)
            for ligand_directory in self.ligand_directories:
                if ligand_directory == "Lf":
                    self.foreign_file_list.append("{0}/{1}/mean_traj".format(dir_name, ligand_directory))
                else:
                    self.self_file_list.append("{0}/{1}/mean_traj".format(dir_name, ligand_directory))

        self.compute_hopfield_error()
        compute_output(self.foreign_file_list)
        compute_output(self.self_file_list)


class ParameterTestingSecondOrder(ParameterTesting):
    def __init__(self, arguments=None):
        ParameterTesting.__init__(self, arguments=arguments)
        self.parameter_list = [float(i) for i in np.sort(lognuniform(low=-4, high=2, size=20))]

        if self.arguments.ss or self.arguments.test:
            x = self.parameter_list[:1]
            y = self.parameter_list[-1:]
            self.parameter_list = [x[0], y[0]]

    def add_steps(self, kp_ligand):
        if self.arguments.steps > 0:
            kp_ligand.ligand.add_cycle_1()
        if self.arguments.steps > 1:
            kp_ligand.ligand.add_cycle_2()
        if self.arguments.steps > 2:
            kp_ligand.ligand.add_cycle_3()
        if self.arguments.steps > 3:
            kp_ligand.ligand.add_cycle_4()

    @staticmethod
    def change_parameter(kp_ligand, parameter):
        kp_ligand.ligand.rate_constants.k_negative_loop = parameter


class ParameterTestingSecondOrderNumbers(ParameterTestingSecondOrder):
    def __init__(self, arguments=None):
        ParameterTestingSecondOrder.__init__(self, arguments=arguments)
        self.parameter_list = [int(i) for i in np.sort(lognuniform(low=2, high=4, size=20))]

        if self.arguments.ss or self.arguments.test:
            x = self.parameter_list[:1]
            y = self.parameter_list[-1:]
            self.parameter_list = [x[0], y[0]]

    @staticmethod
    def change_parameter(kp_ligand, parameter):
        kp_ligand.ligand.n_initial["Zap"] = parameter


def compute_output(file_list):
    output_array = []
    for path in file_list:
        try:
            trajectory = np.loadtxt(path)
            if "parameter_0" in path:
                print(path)
            output_array.append(trajectory[-1])
        except:
            print("Error at {0}".format(path))
    if "Lf" in file_list[0]:
        np.savetxt("Lf_output", output_array, fmt='%f')
    elif "Ls" in file_list[0]:
        np.savetxt("Ls_output", output_array, fmt='%f')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Submitting job for calculating P(C0) as function of steps",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--run', action='store_true', default=False,
                        help='Flag for submitting simulations.')
    parser.add_argument('--p_test', action='store_true', default=False,
                        help='Flag for submitting parameter testing simulations.')
    parser.add_argument('--ss', action='store_true', default=False,
                        help="flag for checking if sims approach steady-state.")
    parser.add_argument('--test', action='store_true', default=False,
                        help="flag for testing.")
    parser.add_argument('--steps', dest='steps', action='store', type=int, default=3,
                        help="number of KP steps.")

    args = parser.parse_args()
    print(str(args.steps))

    directory_name = "{0}_step_check".format(args.steps)
    make_and_cd(directory_name)

    kp_parameter_testing = ParameterTestingSecondOrder(arguments=args)
    # kp_parameter_testing.process_output()

    if args.p_test:
        kp_parameter_testing.parameter_testing(run=args.run)
    else:
        kp_parameter_testing.simple_kp(run=args.run)
