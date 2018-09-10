#!/usr/bin/python
import argparse
import os
import subprocess
import time

import numpy as np

from two_species import KPSingleSpecies


class BindingParameters(object):
    def __init__(self):
        # Zeroth Cycle ligand binding
        self.k_L_on = 0.0022
        self.k_foreign_off = 0.2
        self.k_self_off = 10.0 * self.k_foreign_off

        # First Cycle Lck binding
        self.k_lck_on_R_pmhc = 0.001  # / 200 # self.k_foreign_off / 100.0
        self.k_lck_off_R_pmhc = self.k_foreign_off / 40.0

        self.k_lck_on_R = 1.0 / 100000.0
        self.k_lck_off_R = 20.0

        # Second Cycle Lck phosphorylation
        self.k_p_on_R_pmhc = 2.0
        self.k_p_off_R_pmhc = self.k_foreign_off / 40.0  # self.k_foreign_off / 10.0

        self.k_lck_on_RP = 1.0 / 100000.0
        self.k_lck_off_RP = 20.0

        self.k_p_on_R = 1.0 / 100000.0
        self.k_p_off_R = 20.0

        # Third Cycle Zap binding
        self.k_zap_on_R_pmhc = 0.001  # self.k_lck_on_R_pmhc # self.k_foreign_off / 100.0  # self.k_lck_on_R_pmhc
        self.k_zap_off_R_pmhc = self.k_lck_off_R_pmhc

        self.k_lck_on_zap_R = self.k_lck_on_RP
        self.k_lck_off_zap_R = 0.3

        self.k_zap_on_R = 1.0 / 100000.0
        self.k_zap_off_R = 20.0

        # Fourth Cycle phosphorylate zap
        self.k_p_on_zap_species = self.k_p_on_R_pmhc
        self.k_p_off_zap_species = self.k_p_off_R_pmhc

        self.k_lck_on_zap_p = self.k_lck_on_RP
        self.k_lck_off_zap_p = self.k_lck_off_zap_R

        # Fifth Negative Feedback Loop

        self.k_negative_loop = 20.0

        # # Sixth LAT on
        self.k_lat_on_species = 0.001 / 200  # self.k_lck_on_R_pmhc / 200 # self.k_lck_on_R_pmhc / 100.0  # self.k_foreign_off / 100.0
        self.k_lat_off_species = self.k_lck_off_R_pmhc

        self.k_lat_on_rp_zap = self.k_zap_on_R
        self.k_lat_off_rp_zap = 20.0

        # Seventh LAT phosphorylation
        self.k_p_lat_on_species = 0.008  # self.k_p_on_zap_species / 100
        self.k_p_lat_off_species = self.k_p_off_R_pmhc

        # Eighth LAT -> final product
        self.k_product_on = 0.008
        self.k_product_off = self.k_p_off_R_pmhc

        # Ninth: Positive Feedback Loop
        self.k_positive_loop = 0.0025

        self.k_p_on_lat = self.k_p_on_R_pmhc / 10000.0
        self.k_p_off_lat = 20.0

        # # Eighth Grb on
        # self.k_grb_on_species = self.k_lck_on_R_pmhc
        # self.k_grb_off_species = self.k_lck_off_R_pmhc
        #
        # # Ninth Sos on
        # self.k_sos_on_species = self.k_lck_on_R_pmhc
        # self.k_sos_off_species = self.k_lck_off_R_pmhc
        #
        # # Tenth Ras-GDP on
        # self.k_ras_on_species = self.k_lck_on_R_pmhc
        # self.k_ras_off_species = self.k_lck_off_R_pmhc


class TcrSelfWithForeign(object):
    def __init__(self, arguments=None):
        self.rate_constants = BindingParameters()
        self.arguments = arguments

        # self.n_initial = {"R": 10000, "Lck": 10000, "Zap": 2300000, "Lat": 3000, "Plc": 460000, "Pip": 45000,
        #                   "Lf": 20, "Ls": 0}

        self.n_initial = {"R": 10000, "Lck": 10000, "Zap": 10000, "Lf": 20, "Ls": 0}
        self.record = ["Lf", "Ls", "C_RL", "D_RL"]

        self.simulation_name = "kp_competing"

        self.forward_rates = {"RLf": self.rate_constants.k_L_on, "RLs": self.rate_constants.k_L_on}
        self.reverse_rates = {"C_RL": self.rate_constants.k_foreign_off, "D_RL": self.rate_constants.k_self_off}

        self.forward_rxns = [[["R", "Lf"], ["C_RL"]], [["R", "Ls"], ["D_RL"]]]
        self.reverse_rxns = [[["C_RL"], ["R", "Lf"]], [["D_RL"], ["R", "Ls"]]]

        self.mu = 6
        self.sigma = 1.0
        self.num_samples = 0
        self.set_num_samples()

        self.p_ligand = [int(i) for i in np.round(np.random.lognormal(self.mu, self.sigma, self.num_samples))]

        self.num_kp_steps = 0
        self.output = ["C", "D"]

    def set_num_samples(self):
        if self.arguments:
            if self.arguments.run:
                self.num_samples = 1000
            else:
                self.num_samples = 2

    def change_ligand_concentration(self, concentration):
        self.n_initial["Ls"] = concentration

    def increment_step(self):
        self.num_kp_steps += 1
        self.simulation_name = "kp_steps_" + str(self.num_kp_steps)
        print("Num KP steps = " + str(self.num_kp_steps))


class TcrCycleSelfWithForeign(TcrSelfWithForeign):
    def __init__(self, arguments=None):
        TcrSelfWithForeign.__init__(self, arguments=arguments)

        self.n_initial = {"R": 10000, "Lck": 2000, "Zap": 2000, "LAT": 2000, "Lf": self.arguments.ls_lf, "Ls": 1000}
        # "Grb": 150, "Sos": 150, "Ras_GDP": 100,
        self.record = ["Lf", "Ls"]
        self.output = ["Lf", "Ls"]

        self.simulation_name = "kp_steps_" + str(self.num_kp_steps)

        self.k_L_off = {"Lf": self.rate_constants.k_foreign_off, "Ls": self.rate_constants.k_self_off}

        self.forward_rates = {}
        self.reverse_rates = {}

        self.forward_rxns = []
        self.reverse_rxns = []

    def modify_forward_reverse(self, reactants, products, forward_rate, reverse_rate):
        forward_key = ''.join(reactants) + '_' + ''.join(products)
        self.forward_rates[forward_key] = forward_rate
        self.forward_rxns.append([reactants, products])

        reverse_key = ''.join(products) + '_' + ''.join(reactants)
        self.reverse_rates[reverse_key] = reverse_rate
        self.reverse_rxns.append([products, reactants])

    def create_loop(self, reactants, products, rate):
        forward_key = ''.join(reactants) + '_' + ''.join(products)
        self.forward_rates[forward_key] = rate
        self.forward_rxns.append([reactants, products])

    def add_step_0(self):
        for i in self.output:
            self.modify_forward_reverse(["R", i], ["R" + i], self.rate_constants.k_L_on, self.k_L_off[i])
            self.record.append(''.join(["R" + i]))

    def ligand_off(self, final_product, i):
        new = final_product.replace(i, "")
        self.modify_forward_reverse([final_product], [new, i], self.k_L_off[i],
                                    self.rate_constants.k_L_on)
        return new

    def lck_off(self, no_ligand_product, k_lck_off):
        new = no_ligand_product.replace("_Lck", "")
        self.modify_forward_reverse([no_ligand_product], [new, "Lck"], k_lck_off,
                                    self.rate_constants.k_lck_on_R)
        return new

    def zap_off(self, product):
        new = product.replace("_Zap", "")
        self.modify_forward_reverse([product], [new, "Zap"], self.rate_constants.k_zap_off_R,
                                    self.rate_constants.k_zap_on_R)
        return new

    def dephosphorylate_zap(self, product):
        new = product.replace("_Zap_P", "_Zap")
        self.modify_forward_reverse([product], [new], self.rate_constants.k_p_off_R,
                                    self.rate_constants.k_p_on_R)
        return new

    def rp_zap_off(self, product):
        new = product.replace("RP_Zap_", "")
        self.modify_forward_reverse([product], [new, "RP_Zap"], self.rate_constants.k_lat_off_rp_zap,
                                    self.rate_constants.k_lat_on_rp_zap)
        return new

    def dephosphorylate_lat(self, product):
        new = product.replace("LATP", "LAT")
        self.modify_forward_reverse([product], [new], self.rate_constants.k_p_off_lat,
                                    self.rate_constants.k_p_on_lat)
        return new

    def grb_off(self, product):
        new = product.replace("_Grb", "")
        self.modify_forward_reverse([product], [new, "Grb"], self.rate_constants.k_zap_off_R,
                                    self.rate_constants.k_zap_on_R)
        return new

    def sos_off(self, product):
        new = product.replace("_Sos", "")
        self.modify_forward_reverse([product], [new, "Sos"], self.rate_constants.k_zap_off_R,
                                    self.rate_constants.k_zap_on_R)
        return new

    def record_append(self, final_product):
        self.record.append(''.join([final_product]))

    def cycle_1(self, i, count):
        final_product = "R" + i + "_Lck"
        self.modify_forward_reverse(["R" + i, "Lck"], [final_product], self.rate_constants.k_lck_on_R_pmhc,
                                    self.rate_constants.k_lck_off_R_pmhc)
        no_ligand = self.ligand_off(final_product, i)

        if count == 0:
            self.lck_off(no_ligand, self.rate_constants.k_lck_off_R)

        return final_product

    def cycle_2(self, i, count):
        final_product = "RP" + i + "_Lck"
        self.modify_forward_reverse(["R" + i + "_Lck"], [final_product], self.rate_constants.k_p_on_R_pmhc,
                                    self.rate_constants.k_p_off_R_pmhc)

        no_ligand = self.ligand_off(final_product, i)

        # New pathway back
        if count == 0:
            no_lck = self.lck_off(no_ligand, self.rate_constants.k_lck_off_RP)

            self.modify_forward_reverse([no_lck], ["R"], self.rate_constants.k_p_off_R,
                                        self.rate_constants.k_p_on_R)
        return final_product

    def cycle_3(self, i, count):
        final_product = "RP" + i + "_Lck_Zap"
        self.modify_forward_reverse(["RP" + i + "_Lck", "Zap"], [final_product],
                                    self.rate_constants.k_zap_on_R_pmhc,
                                    self.rate_constants.k_zap_off_R_pmhc)
        no_ligand = self.ligand_off(final_product, i)

        # New pathway back
        if count == 0:
            no_lck = self.lck_off(no_ligand, self.rate_constants.k_lck_off_zap_R)

            self.zap_off(no_lck)

        return final_product

    def cycle_4(self, i, count):
        final_product = "RP" + i + "_Lck_Zap_P"
        self.modify_forward_reverse(["RP" + i + "_Lck_Zap"], [final_product],
                                    self.rate_constants.k_p_on_zap_species, self.rate_constants.k_p_off_zap_species)
        no_ligand = self.ligand_off(final_product, i)

        if count == 0:
            no_lck = self.lck_off(no_ligand, self.rate_constants.k_lck_off_zap_R)
            self.dephosphorylate_zap(no_lck)

        return final_product

    def add_negative_feedback(self):
        self.increment_step()
        count = 0
        for i in self.output:
            self.create_loop(["RP" + i + "_Lck_Zap_P", "R" + i + "_Lck"], ["R" + i + "_LckI", "RP" + i + "_Lck_Zap_P"],
                             self.rate_constants.k_negative_loop)
            self.create_loop(["R" + i + "_LckI"], ["R" + i, "LckI"], 20.0)
            if count == 0:
                self.create_loop(["LckI"], ["Lck"], 20.0)

            count += 1

    def cycle_6(self, i, count):
        final_product = "RP" + i + "_Lck_Zap_P_LAT"
        self.modify_forward_reverse(["RP" + i + "_Lck_Zap_P", "LAT"], [final_product],
                                    self.rate_constants.k_lat_on_species, self.rate_constants.k_lat_off_species)

        no_ligand = self.ligand_off(final_product, i)

        if count == 0:
            no_lck = self.lck_off(no_ligand, self.rate_constants.k_lck_off_zap_R)
            no_zap_p = self.dephosphorylate_zap(no_lck)
            self.rp_zap_off(no_zap_p)

        # self.record_append("RP_Zap_P_LAT")
        # self.record_append("RP_Lck_Zap_P_LAT")

        return final_product

    def add_step_7(self):
        self.increment_step()
        for i in self.output:
            final_product = "LATP"
            self.create_loop(["RP" + i + "_Lck_Zap_P_LAT"], ["RP" + i + "_Lck_Zap_P", final_product],
                             self.rate_constants.k_p_lat_on_species)
            self.create_loop([final_product], ["LAT"], self.rate_constants.k_p_lat_off_species)
            self.record_append(final_product)

    def add_step_8(self):
        self.increment_step()
        final_product = "final_product"
        self.modify_forward_reverse(["LATP"], [final_product], self.rate_constants.k_product_on,
                                    self.rate_constants.k_product_off)
        self.record_append(final_product)

        # self.increment_step()
        # for i in self.output:
        #     final_product = i + "_Product"
        #     self.modify_forward_reverse([i + "_LATP"], [final_product], self.rate_constants.k_product_on,
        #                                 self.rate_constants.k_product_off)
        #     self.record_append(final_product)

    def add_nonspecific_positive_fb(self):
        self.increment_step()
        self.create_loop(["LATP", "final_product"], ["final_product", "final_product"],
                         self.rate_constants.k_positive_loop)

    def add_positive_feedback(self):
        self.increment_step()
        for i in self.output:
            self.create_loop([i + "_LATP", i + "_Product"], [i + "_Product", i + "_Product"],
                             self.rate_constants.k_positive_loop)

    def add_cycle(self, cycle):
        self.increment_step()
        count = 0
        for i in self.output:
            final_product = cycle(i, count)

            self.record_append(final_product)
            count += 1


class TcrCycleForeignLigand(TcrCycleSelfWithForeign):
    def __init__(self, arguments=None):
        TcrCycleSelfWithForeign.__init__(self, arguments)

        del self.n_initial['Ls']
        self.n_initial["Lf"] = 1000
        self.record = ["Lf"]

        self.output = ["Lf"]


class TcrCycleSelfLigand(TcrCycleSelfWithForeign):
    def __init__(self, arguments=None):
        TcrCycleSelfWithForeign.__init__(self, arguments=arguments)

        del self.n_initial['Lf']
        self.n_initial["Ls"] = 1000
        self.record = ["Ls"]

        self.output = ["Ls"]

    # def change_ligand_concentration(self, concentration):
    #     self.n_initial["Ls"] = 384


class KPRealistic(KPSingleSpecies):
    def __init__(self, self_foreign=False, arguments=None):
        KPSingleSpecies.__init__(self, self_foreign=self_foreign, arguments=arguments)

        if self.self_foreign_flag:
            self.ligand = TcrCycleSelfWithForeign(arguments=arguments)
        else:
            self.ligand = TcrCycleSelfLigand(arguments=arguments)

        self.forward_rates = self.ligand.forward_rates
        self.forward_rxns = self.ligand.forward_rxns

        self.reverse_rates = self.ligand.reverse_rates
        self.reverse_rxns = self.ligand.reverse_rxns

        self.record = self.ligand.record

        self.num_files = 1000
        self.run_time = 1000
        self.simulation_time = self.set_simulation_time()

        self.file_list = []

    def set_simulation_time(self):
        simulation_time = self.run_time * (15.0 / 1000)

        return simulation_time

    def wait_for_simulations(self):
        while True:
            list1 = []
            false_list = []

            for file_path in self.file_list:
                list1.append(os.path.isfile(file_path))

                if not os.path.isfile(file_path):
                    false_list.append(file_path)

            print("{0} files remaining".format(len(false_list)))

            if all(list1):
                # All elements are True. Therefore all the files exist. Run %run commands
                break
            else:
                # At least one element is False. Therefore not all the files exist. Run FTP commands again
                time.sleep(5)  # wait 30 seconds before checking again

    def main_script(self, run=False):
        sample = []
        for i in range(self.ligand.num_samples):
            directory = "sample_" + str(i)
            s = self.ligand.p_ligand[i]
            sample.append(s)
            self.ligand.change_ligand_concentration(s)
            simulation_name = self.ligand.simulation_name

            self.file_list.append("{0}/mean_traj".format(directory))

            os.makedirs(directory)
            print("Made " + directory)
            os.chdir(directory)
            print("Changed into directory: " + str(os.getcwd()))

            self.generate(simulation_name, self.set_time_step())
            if run:
                (stdout, stderr) = subprocess.Popen(["qsub {0}".format("qsub.sh")], shell=True, stdout=subprocess.PIPE,
                                                    cwd=os.getcwd()).communicate()
            os.chdir(self.home_directory)

        np.savetxt("Ligand_concentrations", sample, fmt='%f')
        np.savetxt("Ligand_concentrations_sorted", np.sort(sample), fmt='%f')

        # if run:
        #     self.wait_for_simulations()
        #     output = PlotOutputHist()
        #     output.compute_output()


def make_and_cd(directory_name):
    if not os.path.exists(directory_name):
        os.makedirs(directory_name)
        print("Made " + directory_name)
    os.chdir(directory_name)
    print("Changed into directory: " + str(os.getcwd()))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Submitting job for calculating P(C0) as function of steps",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--run', action='store_true', default=False,
                        help='Flag for submitting simulations.')
    parser.add_argument('--ss', action='store_true', default=False,
                        help="flag for checking if sims approach steady-state.")
    parser.add_argument('--steps', dest='steps', action='store', type=int, default=0,
                        help="number of KP steps.")
    parser.add_argument('--ls', action='store_true', default=False,
                        help="flag for submitting Ls calculations.")
    parser.add_argument('--ls_lf', dest='ls_lf', action='store', type=int, default=30,
                        help="number of foreign ligands.")
    args = parser.parse_args()

    directory_name = "{0}_step".format(args.steps)
    make_and_cd(directory_name)

    if args.ls:
        sub_directory = "Ls"
        make_and_cd(sub_directory)
        kp = KPRealistic(arguments=args)

    elif args.ls_lf:
        sub_directory = "Ls_Lf_{0}".format(args.ls_lf)
        make_and_cd(sub_directory)
        kp = KPRealistic(self_foreign=True, arguments=args)
    else:
        raise Exception("Need to specify Ls or Ls_Lf")

    kp.ligand.add_step_0()

    if args.steps > 0:
        kp.ligand.add_cycle(kp.ligand.cycle_1)
    if args.steps > 1:
        kp.ligand.add_cycle(kp.ligand.cycle_2)
    if args.steps > 2:
        kp.ligand.add_cycle(kp.ligand.cycle_3)
    if args.steps > 3:
        kp.ligand.add_cycle(kp.ligand.cycle_4)
    # if args.steps > 4:
    #     kp.ligand.add_negative_feedback()
    if args.steps > 5:
        kp.ligand.add_cycle(kp.ligand.cycle_6)
    if args.steps > 6:
        kp.ligand.add_step_7()
    if args.steps > 7:
        kp.ligand.add_step_8()
    if args.steps > 8:
        kp.ligand.add_nonspecific_positive_fb()

    kp.main_script(run=args.run)
