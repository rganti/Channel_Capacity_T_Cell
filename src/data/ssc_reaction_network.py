import argparse
import os
import subprocess
import time

import numpy as np

from simulation_parameters import InitialConcentrations, DiffusionRates, BindingParameters
from src.data.self_foreign_modules import KPSingleSpecies
from src.general.directory_handling import make_and_cd


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

        self.num_kp_steps = 0
        self.output = ["C", "D"]

    def change_ligand_concentration(self, concentration):
        self.n_initial["Ls"] = concentration

    def increment_step(self):
        self.num_kp_steps += 1
        self.simulation_name = "kp_steps_" + str(self.num_kp_steps)
        print("Num KP steps = " + str(self.num_kp_steps))


class TcrCycleSelfWithForeign(TcrSelfWithForeign):
    def __init__(self, arguments=None):
        TcrSelfWithForeign.__init__(self, arguments=arguments)

        self.initial = InitialConcentrations()
        self.diffusion_constants = DiffusionRates()

        self.n_initial = {"R": self.initial.r_0, "Lf": self.arguments.ls_lf, "Ls": 1000}

        if self.arguments.steps > 0:
            self.n_initial["Lck"] = self.initial.lck_0

        if self.arguments.steps > 2:
            self.n_initial["Zap"] = self.initial.zap_0

        if self.arguments.steps > 4:
            self.n_initial["LAT"] = self.initial.lat_0

        self.record = ["Lf", "Ls"]
        self.output = ["Lf", "Ls"]

        self.record_flag = False

        self.num_kp_steps = self.arguments.steps
        self.simulation_name = "kp_steps_" + str(self.num_kp_steps)

        self.k_L_off = {"Lf": self.rate_constants.k_foreign_off, "Ls": self.rate_constants.k_self_off}

        self.forward_rates = {}
        self.reverse_rates = {}

        self.forward_rxns = []
        self.reverse_rxns = []

        self.diffusion_flag = False
        self.diffusion_rate_dict = {}
        self.diffusion_loc_dict = {}

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

    def update_diffusion(self, species, rate, location="Plasma"):
        self.diffusion_rate_dict[species] = rate
        self.diffusion_loc_dict[species] = location

    def ligand_off(self, final_product, i):
        new = final_product.replace(i, "")
        self.modify_forward_reverse([final_product], [new, i], self.k_L_off[i],
                                    self.rate_constants.k_L_on)
        return new

    def lck_off(self, no_ligand_product, k_lck_off, cycle=3):
        new = no_ligand_product.replace("_Lck", "")
        if cycle < 3:
            self.modify_forward_reverse([no_ligand_product], [new, "Lck"], k_lck_off,
                                        self.rate_constants.k_lck_on_R)
        else:
            self.modify_forward_reverse([no_ligand_product], [new, "Lck"], k_lck_off,
                                        self.rate_constants.k_lck_on_RP)
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

    def add_step_0(self):
        for i in self.output:
            self.modify_forward_reverse(["R", i], ["R" + i], self.rate_constants.k_L_on, self.k_L_off[i])

            if self.arguments.steps == 0:
                self.record.append(''.join(["R" + i]))

        if self.diffusion_flag:
            self.update_diffusion("R", self.diffusion_constants.d_r)

            for i in self.output:
                self.update_diffusion(i, self.diffusion_constants.d_l)
                self.update_diffusion("R" + i, self.diffusion_constants.d_rl)

    def set_record_flag(self, species, step, count):
        if self.arguments.steps == step:
            self.record_flag = True
            if count == 0:
                self.record.append(species)

    def cycle_1(self, i, count):
        final_product = "R" + i + "_Lck"
        self.modify_forward_reverse(["R" + i, "Lck"], [final_product], self.rate_constants.k_lck_on_R_pmhc,
                                    self.rate_constants.k_lck_off_R_pmhc)
        no_ligand = self.ligand_off(final_product, i)

        self.set_record_flag(no_ligand, 1, count)

        if count == 0:
            self.lck_off(no_ligand, self.rate_constants.k_lck_off_R, cycle=1)

        if self.diffusion_flag:
            if count == 0:
                self.update_diffusion("Lck", self.diffusion_constants.d_lck)
                self.update_diffusion(no_ligand, self.diffusion_constants.d_lck)

            self.update_diffusion(final_product, self.diffusion_constants.d_lck)

        return final_product

    def cycle_2(self, i, count):
        final_product = "RP" + i + "_Lck"
        self.modify_forward_reverse(["R" + i + "_Lck"], [final_product], self.rate_constants.k_p_on_R_pmhc,
                                    self.rate_constants.k_p_off_R_pmhc)

        no_ligand = self.ligand_off(final_product, i)

        # New pathway back
        if count == 0:
            no_lck = self.lck_off(no_ligand, self.rate_constants.k_lck_off_RP, cycle=2)

            self.modify_forward_reverse([no_lck], ["R"], self.rate_constants.k_p_off_R,
                                        self.rate_constants.k_p_on_R)

        self.set_record_flag(no_ligand, 2, count)

        if self.diffusion_flag:
            if count == 0:
                self.update_diffusion(no_ligand, self.diffusion_constants.d_lck)
                self.update_diffusion(no_lck, self.diffusion_constants.d_rl)

            self.update_diffusion(final_product, self.diffusion_constants.d_lck)

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

        self.set_record_flag(no_ligand, 3, count)

        if self.diffusion_flag:
            if count == 0:
                self.update_diffusion("Zap", self.diffusion_constants.d_zap, location="Cytosol")
                # self.record_append("Zap in Plasma")
                # self.record_append("Zap in Cytosol")

                self.update_diffusion(no_ligand, self.diffusion_constants.d_lck)
                self.update_diffusion(no_lck, self.diffusion_constants.d_r)

            self.update_diffusion(final_product, self.diffusion_constants.d_lck)

        return final_product

    def add_step_4(self):
        # self.increment_step()
        final_product = "Zap_P"
        for i in self.output:
            self.create_loop(["RP" + i + "_Lck_Zap"], ["RP" + i + "_Lck", final_product],
                             self.rate_constants.k_p_on_zap_species)

        self.create_loop([final_product], ["Zap"], self.rate_constants.k_p_off_zap_species)
        self.record_append(final_product)

    def cycle_4(self, i, count):
        final_product = "RP" + i + "_Lck_Zap_P"
        self.modify_forward_reverse(["RP" + i + "_Lck_Zap"], [final_product],
                                    self.rate_constants.k_p_on_zap_species, self.rate_constants.k_p_off_zap_species)
        no_ligand = self.ligand_off(final_product, i)

        if self.arguments.steps == 4:
            self.record_flag = True

        if count == 0:
            no_lck = self.lck_off(no_ligand, self.rate_constants.k_lck_off_zap_R)
            self.dephosphorylate_zap(no_lck)

        return final_product

    def add_negative_feedback(self):
        # self.increment_step()
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

        if self.arguments.steps == 6:
            self.record_flag = True

        if count == 0:
            no_lck = self.lck_off(no_ligand, self.rate_constants.k_lck_off_zap_R)
            no_zap_p = self.dephosphorylate_zap(no_lck)
            self.rp_zap_off(no_zap_p)

        # self.record_append("RP_Zap_P_LAT")
        # self.record_append("RP_Lck_Zap_P_LAT")

        return final_product

    def cycle_7(self, i, count):
        final_product = "RP{0}_Lck_Zap_P_LATP".format(i)

        self.modify_forward_reverse(["RP" + i + "_Lck_Zap_P_LAT"], [final_product],
                                    self.rate_constants.k_p_lat_1, self.rate_constants.k_p_lat_off_species)

        no_ligand = self.ligand_off(final_product, i)

        if self.arguments.steps == 7:
            self.record_flag = True

        if count == 0:
            no_lck = self.lck_off(no_ligand, self.rate_constants.k_lck_off_zap_R)
            no_zap_p = self.dephosphorylate_zap(no_lck)
            no_latp = self.dephosphorylate_lat(no_zap_p)

        return final_product

    def add_step_8(self):
        final_product = "LATPP"
        for i in self.output:
            self.create_loop(["RP{0}_Lck_Zap_P_LATP".format(i)], ["RP{0}_Lck_Zap_P".format(i), final_product],
                             self.rate_constants.k_p_lat_2)

        self.create_loop([final_product], ["LAT"], self.rate_constants.k_p_lat_off_species)
        self.record_append(final_product)

        return final_product

    def cycle_8(self, i, count):
        final_product = "RP{0}_Lck_Zap_P_LATPP".format(i)

        self.modify_forward_reverse(["RP" + i + "_Lck_Zap_P_LATP"], [final_product],
                                    self.rate_constants.k_p_lat_2, self.rate_constants.k_p_lat_off_species)

        no_ligand = self.ligand_off(final_product, i)

        if self.arguments.steps == 8:
            self.record_flag = True

        if count == 0:
            no_lck = self.lck_off(no_ligand, self.rate_constants.k_lck_off_zap_R)
            no_zap_p = self.dephosphorylate_zap(no_lck)
            no_latp = self.dephosphorylate_lat(no_zap_p)

        return final_product

    # def add_step_7(self):
    #     # self.increment_step()
    #     final_product = "LATP"
    #     for i in self.output:
    #         self.create_loop(["RP" + i + "_Lck_Zap_P_LAT"], ["RP" + i + "_Lck_Zap_P", final_product],
    #                          self.rate_constants.k_p_lat_on_species)
    #
    #     self.create_loop([final_product], ["LAT"], self.rate_constants.k_p_lat_off_species)
    #     self.record_append(final_product)

    # def add_step_8(self):
    #     # self.increment_step()
    #     final_product = "Final_Product"
    #     self.modify_forward_reverse(["LATP"], [final_product], self.rate_constants.k_product_on,
    #                                 self.rate_constants.k_product_off)
    #     self.record_append(final_product)
    #
    #     # self.increment_step()
    #     # for i in self.output:
    #     #     final_product = i + "_Product"
    #     #     self.modify_forward_reverse([i + "_LATP"], [final_product], self.rate_constants.k_product_on,
    #     #                                 self.rate_constants.k_product_off)
    #     #     self.record_append(final_product)

    def add_nonspecific_positive_fb(self):
        # self.increment_step()
        final_product = "Final_Product"
        self.create_loop(["LATP", final_product], [final_product, final_product],
                         self.rate_constants.k_positive_loop)

    # def add_positive_feedback(self):
    #     self.increment_step()
    #     for i in self.output:
    #         self.create_loop([i + "_LATP", i + "_Product"], [i + "_Product", i + "_Product"],
    #                          self.rate_constants.k_positive_loop)

    def add_cycle(self, cycle):
        # self.increment_step()
        count = 0
        for i in self.output:
            final_product = cycle(i, count)

            if self.record_flag:
                self.record_append(final_product)
            count += 1


class ZapDissociation(TcrCycleSelfWithForeign):
    def __init__(self, arguments=None):
        TcrCycleSelfWithForeign.__init__(self, arguments=arguments)

    def add_step_4(self):
        final_product = "Zap_P"
        for i in self.output:
            self.create_loop(["RP{0}_Lck_Zap".format(i)], ["RP{0}_Lck".format(i), final_product],
                             self.rate_constants.k_p_on_zap_species)

        self.create_loop([final_product], ["Zap"], self.rate_constants.k_p_off_zap_species)
        self.record_append(final_product)

    def add_step_6(self):
        final_product = "Zap_P_LAT"
        self.modify_forward_reverse(["Zap_P", "LAT"], [final_product], self.rate_constants.k_lat_on_species,
                                    self.rate_constants.k_lat_off_species)

        self.create_loop([final_product], ["Zap_LAT"], self.rate_constants.k_p_off_R_pmhc)
        self.create_loop(["Zap_LAT"], ["Zap", "LAT"], 20.0)
        self.record_append(final_product)

    def add_step_7(self):
        final_product = "Zap_P_LATP"
        self.modify_forward_reverse(["Zap_P_LAT"], [final_product], self.rate_constants.k_p_lat_1,
                                    self.rate_constants.k_p_lat_off_species)
        self.create_loop([final_product], ["Zap_LATP"], self.rate_constants.k_p_off_zap_species)
        self.create_loop(["Zap_LATP"], ["Zap", "LATP"], 20.0)
        self.create_loop(["LATP"], ["LAT"], self.rate_constants.k_p_lat_off_species)

        self.record_append(final_product)

    def add_step_8(self):
        final_product = "Zap_P_LATPP"
        self.modify_forward_reverse(["Zap_P_LATP"], [final_product], self.rate_constants.k_p_lat_2,
                                    self.rate_constants.k_p_lat_off_species)
        self.create_loop([final_product], ["Zap_LATPP"], self.rate_constants.k_p_off_zap_species)
        self.create_loop(["Zap_LATPP"], ["Zap", "LATPP"], 20.0)
        self.create_loop(["LATPP"], ["LAT"], self.rate_constants.k_p_lat_off_species)

        self.record_append(final_product)


class SelfZapDissociation(ZapDissociation):
    def __init__(self, arguments=None):
        ZapDissociation.__init__(self, arguments=arguments)

        del self.n_initial['Lf']
        self.n_initial["Ls"] = 1000
        self.record = ["Ls"]
        self.output = ["Ls"]


class SscPositiveFbLoop(TcrCycleSelfWithForeign):
    def __init__(self, arguments=None):
        TcrCycleSelfWithForeign.__init__(self, arguments=arguments)
        self.n_initial["Sos"] = 500
        self.n_initial["Ras_GDP"] = 300
        self.n_initial["Ras_GAP"] = 10

    def add_step_8(self):
        # self.increment_step()
        final_product = "LATP_Sos"
        self.modify_forward_reverse(["LATP", "Sos"], [final_product], self.rate_constants.k_sos_on,
                                    self.rate_constants.k_sos_off)
        self.record_append(final_product)

    def add_step_9(self):
        # self.increment_step()
        final_product_1 = "LATP_Sos_Ras_GDP"
        self.modify_forward_reverse(["LATP_Sos", "Ras_GDP"], [final_product_1], self.rate_constants.k_sos_on_rgdp,
                                    self.rate_constants.k_sos_off_rgdp)

        final_product_2 = "LATP_Sos_Ras_GTP"
        self.modify_forward_reverse(["LATP_Sos", "Ras_GTP"], [final_product_2], self.rate_constants.k_sos_on_rgtp,
                                    self.rate_constants.k_sos_off_rgtp)

    def add_step_10(self):
        # self.increment_step()
        previous_product_1 = "LATP_Sos_Ras_GTP"
        final_product_1 = previous_product_1 + "_Ras_GDP"
        self.modify_forward_reverse([previous_product_1, "Ras_GDP"], [final_product_1],
                                    self.rate_constants.k_rgdp_on_sos_rgtp,
                                    self.rate_constants.k_rgdp_off_sos_rgtp)
        self.create_loop([final_product_1], [previous_product_1, "Ras_GTP"], self.rate_constants.k_cat_3)

        previous_product_2 = "LATP_Sos_Ras_GDP"
        final_product_2 = previous_product_2 + "_Ras_GDP"
        self.modify_forward_reverse([previous_product_2, "Ras_GDP"], [final_product_2],
                                    self.rate_constants.k_rgdp_on_sos_rgdp,
                                    self.rate_constants.k_rgdp_off_sos_rgdp)
        self.create_loop([final_product_2], [previous_product_2, "Ras_GTP"], self.rate_constants.k_cat_4)

        product_3 = "Ras_GAP_Ras_GTP"
        self.modify_forward_reverse(["Ras_GAP", "Ras_GTP"], [product_3],
                                    self.rate_constants.k_rgap_on_rgtp,
                                    self.rate_constants.k_rgap_off_rgtp)
        self.create_loop([product_3], ["Ras_GAP", "Ras_GDP"], self.rate_constants.k_cat_5)


class TcrCycleForeignLigand(TcrCycleSelfWithForeign):
    def __init__(self, arguments=None):
        TcrCycleSelfWithForeign.__init__(self, arguments)

        del self.n_initial['Ls']
        self.n_initial["Lf"] = 1000
        self.record = ["Lf"]
        self.output = ["Lf"]


class SelfSscPositiveFbLoop(SscPositiveFbLoop):
    def __init__(self, arguments=None):
        SscPositiveFbLoop.__init__(self, arguments=arguments)

        del self.n_initial['Lf']
        self.n_initial["Ls"] = 1000
        self.record = ["Ls"]
        self.output = ["Ls"]


class TcrCycleSelfLigand(TcrCycleSelfWithForeign):
    def __init__(self, arguments=None):
        TcrCycleSelfWithForeign.__init__(self, arguments=arguments)

        del self.n_initial['Lf']
        self.n_initial["Ls"] = 1000
        self.record = ["Ls"]
        self.output = ["Ls"]


class KPRealistic(KPSingleSpecies):
    def __init__(self, self_foreign=False, arguments=None):
        KPSingleSpecies.__init__(self, self_foreign=self_foreign, arguments=arguments)

        if self.self_foreign_flag:
            # self.ligand = SscPositiveFbLoop(arguments=arguments)
            # self.ligand = TcrCycleSelfWithForeign(arguments=arguments)
            self.ligand = ZapDissociation(arguments=arguments)
        else:
            # self.ligand = SelfSscPositiveFbLoop(arguments=arguments)
            # self.ligand = TcrCycleSelfLigand(arguments=arguments)
            self.ligand = SelfZapDissociation(arguments=arguments)

        self.num_files = 200
        self.run_time = 500

        self.mu = 6
        self.sigma = 1.0
        self.num_samples = self.set_num_samples()

        self.p_ligand = [int(i) for i in np.round(np.random.lognormal(self.mu, self.sigma, self.num_samples))]

        self.file_list = []

    def set_num_samples(self):
        if self.arguments.run:
            num_samples = 1000
        else:
            num_samples = 2

        return num_samples

    def set_simulation_time(self, ls=500):
        if ls < 500:
            simulation_time = 4.0
        else:
            simulation_time = 40.0

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

    def create_steps(self):
        self.ligand.add_step_0()

        if self.arguments.steps > 0:
            self.ligand.add_cycle(self.ligand.cycle_1)
        if self.arguments.steps > 1:
            self.ligand.add_cycle(self.ligand.cycle_2)
        if self.arguments.steps > 2:
            self.ligand.add_cycle(self.ligand.cycle_3)
        if self.arguments.steps > 3:
            self.ligand.add_step_4()
        if self.arguments.steps > 5:
            self.ligand.add_step_6()
        if self.arguments.steps > 6:
            self.ligand.add_step_7()
        if self.arguments.steps > 7:
            self.ligand.add_step_8()
        # if self.arguments.steps > 3:
        #     self.ligand.add_cycle(self.ligand.cycle_4)
        # if self.arguments.steps > 5:
        #     self.ligand.add_cycle(self.ligand.cycle_6)
        # if self.arguments.steps > 6:
        #     self.ligand.add_cycle(self.ligand.cycle_7)
        # if self.arguments.steps > 7:
        #     # self.ligand.add_cycle(self.ligand.cycle_8)
        #     self.ligand.add_step_8()

    def main_script(self, run=False):
        sample = []
        self.create_steps()

        for i in range(len(self.p_ligand)):
            time.sleep(0.5)
            directory = "sample_" + str(i)

            s = self.p_ligand[i]
            sample.append(s)
            self.ligand.change_ligand_concentration(s)
            simulation_name = self.ligand.simulation_name

            self.file_list.append("{0}/mean_traj".format(directory))

            os.makedirs(directory)
            print("Made " + directory)
            os.chdir(directory)
            print("Changed into directory: " + str(os.getcwd()))

            self.generate(simulation_name, self.set_time_step(), ls=s)
            if run:
                (stdout, stderr) = subprocess.Popen(["qsub {0}".format("qsub.sh")], shell=True, stdout=subprocess.PIPE,
                                                    cwd=os.getcwd()).communicate()
            os.chdir(self.home_directory)

        print(str(self.ligand.record))
        np.savetxt("Ligand_concentrations", sample, fmt='%f')
        np.savetxt("Ligand_concentrations_sorted", np.sort(sample), fmt='%f')

        # if run:
        #     self.wait_for_simulations()
        #     output = PlotOutputHist()
        #     output.compute_output()


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

    kp.main_script(run=args.run)
