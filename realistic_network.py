#!/usr/bin/python
import argparse
import os

import numpy as np

from two_species import KPSingleSpecies


class BindingParameters(object):
    def __init__(self):
        # Zeroth Cycle ligand binding
        self.k_L_on = 0.0022
        self.k_foreign_off = 0.2
        self.k_self_off = 10.0 * self.k_foreign_off

        # First Cycle Lck binding

        self.k_lck_on_R_pmhc = 10.0
        self.k_lck_off_R_pmhc = self.k_foreign_off / 10.0

        self.k_lck_on_R = self.k_lck_on_R_pmhc / 10000.0
        self.k_lck_off_R = 20.0

        # Second Cycle Lck phosphorylation
        self.k_p_on_R_pmhc = 5.0
        self.k_p_off_R_pmhc = self.k_foreign_off / 10.0

        self.k_lck_on_RP = self.k_lck_on_R_pmhc / 10000.0
        self.k_lck_off_RP = 20.0

        self.k_p_on_R = self.k_p_on_R_pmhc / 1000
        self.k_p_off_R = 10.0

        # all of the above seems to work for 1st order kinetics

        # Third Cycle Zap binding
        self.k_zap_on_R_pmhc = 3 * self.k_lck_on_R_pmhc
        self.k_zap_off_R_pmhc = self.k_foreign_off / 10.0

        self.k_lck_on_zap_R = self.k_lck_on_RP
        self.k_lck_off_zap_R = self.k_L_on / 10.0
        self.k_lck_off_zap_R = 1.0

        self.k_zap_on_R = self.k_zap_on_R_pmhc / 10000.0
        self.k_zap_off_R = 20.0

        # self.k_tcr_tcrp_on = self.k_L_on / 10000.0
        # self.k_tcr_tcrp_off = self.k_p_off_R * 10.0

        # self.k_lat_on = 0.1
        # self.k_lat_off = 0.1
        #
        # self.k_plc_on = 161.2
        # self.k_plc_off = 1.0
        #
        # self.k_p_plc_on = 0.01
        #
        # self.k_pip_on = 0.1
        # self.k_pip_off = 0.2


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
            if self.arguments.test or self.arguments.ss:
                self.num_samples = 5
            else:
                self.num_samples = 1000
        else:
            self.num_samples = 1000

    def change_ligand_concentration(self, concentration):
        self.n_initial["Ls"] = concentration

    def increment_step(self):
        self.num_kp_steps += 1
        self.simulation_name = "kp_steps_" + str(self.num_kp_steps)
        print("Num KP steps = " + str(self.num_kp_steps))


class TcrStepsSelfWithForeign(TcrSelfWithForeign):
    def __init__(self):
        TcrSelfWithForeign.__init__(self)

    def modify_forward_reverse(self, reactants, products, forward_rate, reverse_rate):
        self.forward_rates[''.join(reactants)] = forward_rate
        self.forward_rxns.append([reactants, products])

        self.reverse_rates[''.join(products)] = reverse_rate
        self.reverse_rxns.append([products, reactants])

        self.record.append(''.join(products))

    def add_step_1(self):
        self.increment_step()
        for i in self.output:
            self.modify_forward_reverse([i + "_RL", "Lck"], [i + "_RL1"], self.rate_constants.k_lck_on_R_pmhc,
                                        self.rate_constants.k_lck_off_R_pmhc)

    def add_step_2(self):
        self.increment_step()
        for i in self.output:
            self.modify_forward_reverse([i + "_RL1"], [i + "_RL2"], self.rate_constants.k_p_on,
                                        self.rate_constants.k_p_off)

    def add_step_3(self):
        self.increment_step()
        for i in self.output:
            self.modify_forward_reverse([i + "_RL2", "Zap"], [i + "_Zap"], self.rate_constants.k_zap_on_R_pmhc,
                                        self.rate_constants.k_zap_off_R_pmhc)

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


class TcrSelfLigand(TcrStepsSelfWithForeign):
    def __init__(self):
        TcrStepsSelfWithForeign.__init__(self)

        del self.n_initial['Lf']
        self.record = ["Ls", "D_RL"]
        self.simulation_name = "kp_ls"

        self.forward_rates = {"RLs": self.rate_constants.k_L_on}
        self.reverse_rates = {"D_RL": self.rate_constants.k_self_off}

        self.forward_rxns = [[["R", "Ls"], ["D_RL"]]]
        self.reverse_rxns = [[["D_RL"], ["R", "Ls"]]]

        self.output = ["D"]


class TcrCycleSelfWithForeign(TcrSelfWithForeign):
    def __init__(self, arguments=None):
        TcrSelfWithForeign.__init__(self, arguments=arguments)

        self.n_initial = {"R": 10000, "Lck": 10000, "Zap": 10000, "Lf": 100, "Ls": 0}
        self.record = ["Lf", "Ls"]
        self.output = ["Lf", "Ls"]

        self.simulation_name = "kp_competing"

        self.k_L_off = {"Lf": self.rate_constants.k_foreign_off, "Ls": self.rate_constants.k_self_off}

        self.forward_rates = {}
        self.reverse_rates = {}

        self.forward_rxns = []
        self.reverse_rxns = []

        # # Only foreign
        # del self.n_initial['Ls']
        # self.record = ["Lf"]
        # self.output = ["Lf"]

    def modify_forward_reverse(self, reactants, products, forward_rate, reverse_rate):
        forward_key = ''.join(reactants) + '_' + ''.join(products)
        self.forward_rates[forward_key] = forward_rate
        self.forward_rxns.append([reactants, products])

        reverse_key = ''.join(products) + '_' + ''.join(reactants)
        self.reverse_rates[reverse_key] = reverse_rate
        self.reverse_rxns.append([products, reactants])

    def add_step_0(self):
        # del self.n_initial["Lck"]
        # del self.n_initial["Zap"]
        for i in self.output:
            self.modify_forward_reverse(["R", i], ["R" + i], self.rate_constants.k_L_on,
                                        self.k_L_off[i])
            self.record.append(''.join(["R" + i]))

    def add_cycle_1(self):
        self.increment_step()
        for i in self.output:
            self.modify_forward_reverse(["R" + i, "Lck"], ["R" + i + "_Lck"], self.rate_constants.k_lck_on_R_pmhc,
                                        self.rate_constants.k_lck_off_R_pmhc)

            self.modify_forward_reverse(["R" + i + "_Lck"], ["R_Lck", i], self.k_L_off[i],
                                        self.rate_constants.k_L_on)

            self.modify_forward_reverse(["R_Lck"], ["R", "Lck"], self.rate_constants.k_lck_off_R,
                                        self.rate_constants.k_lck_on_R)

            # self.record.append(''.join(["R_Lck"]))
            self.record.append(''.join(["R" + i + "_Lck"]))

    # def add_cycle_1_first_order(self):
    #     self.increment_step()
    #     for i in self.output:
    #         self.modify_forward_reverse(["R" + i], ["R" + i + "_Lck"], self.rate_constants.k_lck_on_R_pmhc,
    #                                     self.rate_constants.k_lck_off_R_pmhc)
    #
    #         self.modify_forward_reverse(["R" + i + "_Lck"], ["R_Lck", i], self.k_L_off[i],
    #                                     self.rate_constants.k_L_on)
    #
    #         self.modify_forward_reverse(["R_Lck"], ["R"], self.rate_constants.k_lck_off_R,
    #                                     self.rate_constants.k_lck_on_R)
    #
    #         # self.record.append(''.join(["R_Lck"]))
    #         self.record.append(''.join(["R" + i + "_Lck"]))
    #
    # def add_cycle_2_first_order(self):
    #     self.increment_step()
    #     for i in self.output:
    #         self.modify_forward_reverse(["R" + i + "_Lck"], ["RP" + i + "_Lck"], self.rate_constants.k_p_on_R_pmhc,
    #                                     self.rate_constants.k_p_off_R_pmhc)
    #
    #         self.modify_forward_reverse(["RP" + i + "_Lck"], ["RP_Lck", i], self.k_L_off[i],
    #                                     self.rate_constants.k_L_on)
    #
    #         # New pathway back
    #         self.modify_forward_reverse(["RP_Lck"], ["RP"], self.rate_constants.k_lck_off_RP,
    #                                     self.rate_constants.k_lck_on_RP)
    #
    #         self.modify_forward_reverse(["RP"], ["R"], self.rate_constants.k_p_off_R,
    #                                     self.rate_constants.k_p_on_R)
    #
    #         # self.record.append(''.join(["RP"]))
    #         # self.record.append(''.join(["RP_Lck"]))
    #         self.record.append(''.join(["RP" + i + "_Lck"]))

    def add_cycle_2(self):
        self.increment_step()
        for i in self.output:
            self.modify_forward_reverse(["R" + i + "_Lck"], ["RP" + i + "_Lck"], self.rate_constants.k_p_on_R_pmhc,
                                        self.rate_constants.k_p_off_R_pmhc)

            self.modify_forward_reverse(["RP" + i + "_Lck"], ["RP_Lck", i], self.k_L_off[i],
                                        self.rate_constants.k_L_on)

            # New pathway back
            self.modify_forward_reverse(["RP_Lck"], ["RP", "Lck"], self.rate_constants.k_lck_off_RP,
                                        self.rate_constants.k_lck_on_RP)

            self.modify_forward_reverse(["RP"], ["R"], self.rate_constants.k_p_off_R,
                                        self.rate_constants.k_p_on_R)

            # self.modify_forward_reverse(["RP_Lck"], ["R_Lck"], self.rate_constants.k_p_off_R,
            #                             self.rate_constants.k_p_on_R)
            #
            # self.record.append(''.join(["RP"]))
            # self.record.append(''.join(["RP_Lck"]))
            self.record.append(''.join(["RP" + i + "_Lck"]))

    def add_cycle_3(self):
        self.increment_step()
        for i in self.output:
            self.modify_forward_reverse(["RP" + i + "_Lck", "Zap"], ["RP" + i + "_Lck_Zap"],
                                        self.rate_constants.k_zap_on_R_pmhc,
                                        self.rate_constants.k_zap_off_R_pmhc)
            self.modify_forward_reverse(["RP" + i + "_Lck_Zap"], ["RP_Lck_Zap", i],
                                        self.k_L_off[i], self.rate_constants.k_L_on)

            # New pathway back
            self.modify_forward_reverse(["RP_Lck_Zap"], ["RP_Zap", "Lck"], self.rate_constants.k_lck_off_zap_R,
                                        self.rate_constants.k_lck_on_zap_R)

            self.modify_forward_reverse(["RP_Zap"], ["RP", "Zap"], self.rate_constants.k_zap_off_R,
                                        self.rate_constants.k_zap_on_R)

            self.modify_forward_reverse(["RP"], ["R"], self.rate_constants.k_p_off_R,
                                        self.rate_constants.k_p_on_R)

            # self.record.append(''.join(["RP_Lck_Zap"]))
            self.record.append(''.join(["RP" + i + "_Lck_Zap"]))

    # def add_new_cycle_3(self):
    #     self.increment_step()
    #     for i in self.output:
    #         self.modify_forward_reverse(["RP" + i + "_Lck", "Zap"], ["RP" + i + "_Lck_Zap"],
    #                                     self.rate_constants.k_zap_on_R_pmhc,
    #                                     self.rate_constants.k_zap_off_R_pmhc)
    #         self.modify_forward_reverse(["RP" + i + "_Lck_Zap"], ["RP_Lck_Zap", i],
    #                                     self.k_L_off[i], self.rate_constants.k_L_on)
    #
    #         # New pathway back
    #         self.modify_forward_reverse(["RP_Lck_Zap"], ["RP_Zap", "Lck"], self.rate_constants.k_lck_off_zap_R,
    #                                     self.rate_constants.k_lck_on_zap_R)
    #
    #         self.modify_forward_reverse(["RP_Zap"], ["R_Zap"], self.rate_constants.k_p_off_R,
    #                                     self.rate_constants.k_p_on_R)
    #
    #         self.modify_forward_reverse(["R_Zap"], ["R", "Zap"], self.rate_constants.k_zap_off_R,
    #                                     self.rate_constants.k_zap_on_R)
    #
    #         self.record.append(''.join(["RP_Lck_Zap"]))
    #         self.record.append(''.join(["RP" + i + "_Lck_Zap"]))


class TcrCycleSelfWithForeignFirstOrder(TcrCycleSelfWithForeign):
    def __init__(self, arguments=None):
        TcrCycleSelfWithForeign.__init__(self, arguments)
        self.n_initial = {"R": 10000, "Lf": 100, "Ls": 0}

    def add_cycle_1_first_order(self):
        self.increment_step()
        for i in self.output:
            self.modify_forward_reverse(["R" + i], ["R" + i + "_Lck"], self.rate_constants.k_lck_on_R_pmhc,
                                        self.rate_constants.k_lck_off_R_pmhc)

            self.modify_forward_reverse(["R" + i + "_Lck"], ["R_Lck", i], self.k_L_off[i],
                                        self.rate_constants.k_L_on)

            self.modify_forward_reverse(["R_Lck"], ["R"], self.rate_constants.k_lck_off_R,
                                        self.rate_constants.k_lck_on_R)

            self.record.append(''.join(["R" + i + "_Lck"]))

    def add_cycle_2_first_order(self):
        self.increment_step()
        for i in self.output:
            self.modify_forward_reverse(["R" + i + "_Lck"], ["RP" + i + "_Lck"], self.rate_constants.k_p_on_R_pmhc,
                                        self.rate_constants.k_p_off_R_pmhc)

            self.modify_forward_reverse(["RP" + i + "_Lck"], ["RP_Lck", i], self.k_L_off[i],
                                        self.rate_constants.k_L_on)

            # New pathway back
            self.modify_forward_reverse(["RP_Lck"], ["RP"], self.rate_constants.k_lck_off_RP,
                                        self.rate_constants.k_lck_on_RP)

            self.modify_forward_reverse(["RP"], ["R"], self.rate_constants.k_p_off_R,
                                        self.rate_constants.k_p_on_R)

            self.record.append(''.join(["RP" + i + "_Lck"]))

    def add_cycle_3_first_order(self):
        self.increment_step()
        for i in self.output:
            self.modify_forward_reverse(["RP" + i + "_Lck"], ["RP" + i + "_Lck_Zap"],
                                        self.rate_constants.k_zap_on_R_pmhc,
                                        self.rate_constants.k_zap_off_R_pmhc)
            self.modify_forward_reverse(["RP" + i + "_Lck_Zap"], ["RP_Lck_Zap", i],
                                        self.k_L_off[i], self.rate_constants.k_L_on)

            # New pathway back
            self.modify_forward_reverse(["RP_Lck_Zap"], ["RP_Zap"], self.rate_constants.k_lck_off_zap_R,
                                        self.rate_constants.k_lck_on_zap_R)

            self.modify_forward_reverse(["RP_Zap"], ["RP"], self.rate_constants.k_zap_off_R,
                                        self.rate_constants.k_zap_on_R)

            self.record.append(''.join(["RP" + i + "_Lck_Zap"]))


class TcrCycleForeignLigand(TcrCycleSelfWithForeignFirstOrder):
    def __init__(self, arguments=None):
        TcrCycleSelfWithForeignFirstOrder.__init__(self, arguments)

        del self.n_initial['Ls']
        self.n_initial["Lf"] = 1000
        self.record = ["Lf"]
        self.simulation_name = "kp_lf"

        self.output = ["Lf"]


class TcrCycleSelfLigand(TcrCycleSelfWithForeignFirstOrder):
    def __init__(self, arguments=None):
        TcrCycleSelfWithForeignFirstOrder.__init__(self, arguments=arguments)

        del self.n_initial['Lf']
        self.n_initial["Ls"] = 1000
        self.record = ["Ls"]
        self.simulation_name = "kp_ls"

        self.output = ["Ls"]

    # def change_ligand_concentration(self, concentration):
    #     self.n_initial["Ls"] = 264


class KPRealistic(KPSingleSpecies):
    def __init__(self, self_foreign=False, arguments=None):
        KPSingleSpecies.__init__(self, self_foreign=self_foreign, arguments=arguments)

        if self.self_foreign_flag:
            self.ligand = TcrCycleSelfWithForeignFirstOrder(arguments=arguments)
        else:
            self.ligand = TcrCycleSelfLigand(arguments=arguments)

        self.forward_rates = self.ligand.forward_rates
        self.forward_rxns = self.ligand.forward_rxns

        self.reverse_rates = self.ligand.reverse_rates
        self.reverse_rxns = self.ligand.reverse_rxns

        self.record = self.ligand.record

        self.simulation_time = 5


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
        kp = KPRealistic(self_foreign=True, arguments=args)
    elif "Ls" in os.getcwd():
        kp = KPRealistic(arguments=args)
    else:
        raise Exception("Incorrect Directory labeling. Specify (Ls, Ls_Lf)")

    kp.ligand.add_step_0()
    if "1_step" in os.getcwd():
        kp.ligand.add_cycle_1_first_order()
    elif "2_step" in os.getcwd():
        kp.ligand.add_cycle_1_first_order()
        kp.ligand.add_cycle_2_first_order()
    elif "3_step" in os.getcwd():
        kp.ligand.add_cycle_1_first_order()
        kp.ligand.add_cycle_2_first_order()
        kp.ligand.add_cycle_3_first_order()


    print(kp.ligand.forward_rates)
    print(kp.ligand.reverse_rates)

    kp.main_script(run=args.run)
