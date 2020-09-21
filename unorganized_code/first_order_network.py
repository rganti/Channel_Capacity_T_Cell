from realistic_network import TcrSelfWithForeign, TcrCycleSelfWithForeign
from simulation_parameters import BindingParameters


class BindingParametersFirstOrder(BindingParameters):
    def __init__(self):
        BindingParameters.__init__(self)
        self.k_lck_on_R_pmhc = 10.0
        self.k_lck_off_R_pmhc = self.k_foreign_off / 10.0

        self.k_p_off_R_pmhc = self.k_foreign_off / 10.0

        self.k_zap_off_R_pmhc = self.k_foreign_off / 10.0

        self.k_zap_on_R_pmhc = 3 * self.k_lck_on_R_pmhc
        self.k_zap_off_R_pmhc = self.k_foreign_off / 10.0

        self.k_lck_off_zap_R = 1.0


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


class TcrCycleSelfWithForeignFirstOrder(TcrCycleSelfWithForeign):
    def __init__(self, arguments=None):
        TcrCycleSelfWithForeign.__init__(self, arguments)
        self.rate_constants = BindingParametersFirstOrder()
        self.n_initial = {"R": 10000, "Lf": 20, "Ls": 0}

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

    def add_negative_feedback(self):
        self.increment_step()
        for i in self.output:
            self.create_loop(["RP" + i + "_Lck_Zap", "R" + i + "_Lck"], ["R" + i, "RP" + i + "_Lck_Zap"],
                             self.rate_constants.k_negative_loop)


class ExtendedNetwork(TcrCycleSelfWithForeign):
    def __init__(self):
        TcrCycleSelfWithForeign.__init__(self)

    def cycle_7(self, i, count):
        final_product = "RP" + i + "_Lck_Zap_P_LATP"
        self.modify_forward_reverse(["RP" + i + "_Lck_Zap_P_LAT"], [final_product],
                                    self.rate_constants.k_p_lat_on_species, self.rate_constants.k_p_lat_off_species)
        no_ligand = self.ligand_off(final_product, i)

        if count == 0:
            no_lck = self.lck_off(no_ligand, self.rate_constants.k_lck_off_zap_R)
            no_zap_p = self.dephosphorylate_zap(no_lck)
            no_rp_zap = self.rp_zap_off(no_zap_p)
            self.dephosphorylate_lat(no_rp_zap)

        return final_product

    def cycle_8(self, i, count):
        final_product = "RP" + i + "_Lck_Zap_P_LATP_Grb"
        self.modify_forward_reverse(["RP" + i + "_Lck_Zap_P_LATP", "Grb"], [final_product],
                                    self.rate_constants.k_grb_on_species, self.rate_constants.k_grb_off_species)
        no_ligand = self.ligand_off(final_product, i)

        if count == 0:
            no_lck = self.lck_off(no_ligand, self.rate_constants.k_lck_off_zap_R)
            no_zap_p = self.dephosphorylate_zap(no_lck)
            no_rp_zap = self.rp_zap_off(no_zap_p)
            self.grb_off(no_rp_zap)

        return final_product

    def cycle_9(self, i, count):
        final_product = "RP" + i + "_Lck_Zap_P_LATP_Grb_Sos"
        self.modify_forward_reverse(["RP" + i + "_Lck_Zap_P_LATP_Grb", "Sos"], [final_product],
                                    self.rate_constants.k_sos_on_species, self.rate_constants.k_sos_off_species)
        no_ligand = self.ligand_off(final_product, i)

        if count == 0:
            no_lck = self.lck_off(no_ligand, self.rate_constants.k_lck_off_zap_R)
            no_zap_p = self.dephosphorylate_zap(no_lck)
            no_rp_zap = self.rp_zap_off(no_zap_p)
            self.sos_off(no_rp_zap)

        return final_product

    def add_positive_feedback(self):
        self.increment_step()
        count = 0
        for i in self.output:
            # catalytic binding
            self.modify_forward_reverse(["RP" + i + "_Lck_Zap_P_LATP_Grb_Sos", "Ras_GDP"],
                                        ["R" + i + "_species_Sos_Ras_GDP"],
                                        self.rate_constants.k_ras_on_species, self.rate_constants.k_ras_off_species)

            # proofreading of Sos-Ras_GDP complex
            no_ligand = self.ligand_off("R" + i + "_species_Sos_Ras_GDP", i)
            self.modify_forward_reverse([no_ligand], ["RP_Lck_Zap_P_LATP_Grb", "Sos"], self.rate_constants.k_zap_off_R,
                                        self.rate_constants.k_zap_on_R)

            #
            # self.modify_forward_reverse(["RP" + i + "_Lck_Zap_P_LATP_Grb_Sos", "Ras_GDP"],
            #                             ["R" + i + "_species_Sos_Ras_GDP"],
            #                             0.00027, 0.004)

            self.create_loop(["R" + i + "_species_Sos_Ras_GDP"], ["RP" + i + "_Lck_Zap_P_LATP_Grb_Sos", "Ras_GTP"],
                             0.0005)

            # # allosteric binding
            self.modify_forward_reverse(["RP" + i + "_Lck_Zap_P_LATP_Grb_Sos", "Ras_GTP"],
                                        ["R" + i + "_species_Sos_Ras_GTP"], self.rate_constants.k_ras_on_species,
                                        self.rate_constants.k_ras_off_species)

            self.modify_forward_reverse(["R" + i + "_species_Sos_Ras_GTP", "Ras_GDP"],
                                        ["R" + i + "_species_Sos_Ras_GTP_Ras_GDP"],
                                        0.05, 0.1)

            no_ligand = self.ligand_off("R" + i + "_species_Sos_Ras_GTP_Ras_GDP", i)
            self.modify_forward_reverse([no_ligand], ["RP_Lck_Zap_P_LATP_Grb", "Sos"], self.rate_constants.k_zap_off_R,
                                        self.rate_constants.k_zap_on_R)

            self.create_loop(["R" + i + "_species_Sos_Ras_GTP_Ras_GDP"], ["R" + i + "_species_Sos_Ras_GTP", "Ras_GTP"],
                             0.038)

            self.record_append("Ras_GTP")
            count += 1
