'''This script is the main one for generating all possible reactions of main T cell network using the PYSB package.'''

import argparse
import os
import pickle

import numpy as np
from pysb import *
from pysb.integrate import odesolve

from simulation_parameters import InitialConcentrations, BindingParameters


def e(i, s=""):
    eval('{0}{1}()'.format(i, s))


def add_new_monomer(product):
    try:
        model.monomers[product]
    except:
        Monomer(product)
        Initial(eval(product + "()"), Parameter(product + "_0", 0))


def add_observable(species):
    Observable("O_{0}".format(species), eval('{0}()'.format(species)))


def write_columns(observables):
    f = open("column_names", "w")
    for item in observables:
        f.write("{0} ".format(item))
    f.write("\n")
    f.close()


def write_model_attributes(attributes, file_name):
    f = open(file_name, "w")
    for item in attributes:
        f.write("{0} \n".format(item))
    f.write("\n")
    f.close()


class PysbTcrSelfWithForeign(object):
    def __init__(self, steps=8, self_foreign=False, lf=30):
        self.rate_constants = BindingParameters()
        self.initial_conditions = InitialConcentrations()
        self.self_foreign = self_foreign

        self.lf = lf
        self.ls = 1000
        self.steps = steps

        self.run_time = 1000
        self.tspan = np.linspace(0, self.run_time)

        self.mu = 6
        self.sigma = 1.0
        self.num_samples = 1000

        self.p_flag = os.path.exists("parameters.pickle")

        if self.p_flag:
            self.parameters = pickle.load(open("parameters.pickle", "rb"))
            print(self.parameters)

        self.p_ligand = [int(i) for i in np.round(np.random.lognormal(self.mu, self.sigma, self.num_samples))]

        if self.self_foreign:
            self.output = ["Ls", "Lf"]
        else:
            self.output = ["Ls"]

        self.model = Model()

    def define_monomers(self):
        Monomer('R')
        Monomer('Ls')
        Monomer('Lck')
        Monomer('Zap')
        Monomer('LAT')
        Monomer('Grb')
        Monomer('Sos')

        Parameter('R_0', self.initial_conditions.r_0)
        Parameter('Ls_0', self.ls)
        Parameter('Lck_0', self.initial_conditions.lck_0)
        Parameter('Zap_0', self.initial_conditions.zap_0)
        Parameter('LAT_0', self.initial_conditions.lat_0)
        Parameter('Sos_0', self.initial_conditions.sos_0)

        Initial(R(), R_0)
        Initial(Ls(), Ls_0)
        Initial(Lck(), Lck_0)
        Initial(Zap(), Zap_0)
        Initial(LAT(), LAT_0)
        Initial(Sos(), Sos_0)

        # Positive Feedback loop
        Monomer('Ras_GDP')
        Monomer('Ras_GTP')
        Monomer('Ras_GAP')

        Parameter('Ras_GDP_0', self.initial_conditions.ras_gdp_0)
        Parameter('Ras_GAP_0', self.initial_conditions.ras_gap_0)
        Parameter('Ras_GTP_0', 0)

        Initial(Ras_GDP(), Ras_GDP_0)
        Initial(Ras_GAP(), Ras_GAP_0)
        Initial(Ras_GTP(), Ras_GTP_0)

        if "Lf" in self.output:
            Monomer('Lf')
            Parameter('Lf_0', self.lf)
            Initial(Lf(), Lf_0)

    def add_new_monomer(self, product):
        try:
            self.model.monomers[product]
        except:
            Monomer(product)
            Initial(eval(product + "()"), Parameter(product + "_0", 0))

    def define_rates(self):
        Parameter('k_L_on', self.rate_constants.k_L_on)
        if "Lf" in self.output:
            Parameter('k_Lf_off', self.rate_constants.k_foreign_off)
        Parameter('k_Ls_off', self.rate_constants.k_self_off)

    def add_step_0(self, i):
        product = "R{0}".format(i)
        self.add_new_monomer(product)

        Rule('{0}_bind'.format(product), R() + eval('{0}()'.format(i)) | eval('{0}()'.format(product)),
             k_L_on, eval('k_{0}_off'.format(i)))

        add_observable(i)
        add_observable(product)
        return product

    def zap_off(self, product):
        new = product.replace("_Zap", "")
        self.add_new_monomer(new)
        Rule('{0}_unbindZap'.format(product), eval('{0}()'.format(product)) | eval('{0}()'.format(new)) + Zap(),
             k_zap_off_R, k_zap_on_R)
        return new

    def lck_off(self, no_ligand_product, k_off="k_lck_off_R"):
        new = no_ligand_product.replace("_Lck", "")
        self.add_new_monomer(new)
        Rule('{0}_unbindLck'.format(no_ligand_product),
             eval('{0}()'.format(no_ligand_product)) | eval('{0}()'.format(new)) + Lck(),
             eval(k_off), k_lck_on_R)

        return new

    def unbind_ligand(self, i, final_product):
        new = final_product.replace(i, "")
        self.add_new_monomer(new)
        Rule('{0}_unbindL'.format(final_product),
             eval('{0}()'.format(final_product)) | eval('{0}()'.format(i)) + eval('{0}()'.format(new)),
             eval('k_{0}_off'.format(i)), k_L_on)

        return new

    def cycle_1(self, i):
        if i == "Ls":
            if self.p_flag:
                Parameter('k_lck_on_RL', self.parameters['k_lck_on_RL'])
            else:
                Parameter('k_lck_on_RL', self.rate_constants.k_lck_on_R_pmhc)
            Parameter('k_lck_off_RL', self.rate_constants.k_lck_off_R_pmhc)

            Parameter('k_lck_on_R', self.rate_constants.k_lck_on_R)
            Parameter('k_lck_off_R', self.rate_constants.k_lck_off_R)

        previous_product = self.add_step_0(i)

        product = "R{0}_Lck".format(i)
        add_new_monomer(product)
        Rule('{0}_bind'.format(product), eval('{0}()'.format(previous_product)) + Lck() | eval('{0}()'.format(product)),
             k_lck_on_RL, k_lck_off_RL)
        add_observable(product)

        no_ligand = self.unbind_ligand(i, product)
        if i == "Ls":
            lck_off = self.lck_off(no_ligand)

        return product

    def cycle_2(self, i):
        if i == "Ls":
            if self.p_flag:
                Parameter('k_p_on_R_pmhc', self.parameters['k_p_on_R_pmhc'])
            else:
                Parameter('k_p_on_R_pmhc', self.rate_constants.k_p_on_R_pmhc)
            Parameter('k_p_off_R_pmhc', self.rate_constants.k_p_off_R_pmhc)

            Parameter('k_p_on_R', self.rate_constants.k_p_on_R)
            Parameter('k_p_off_R', self.rate_constants.k_p_off_R)

        previous_product = self.cycle_1(i)

        product = "RP{0}_Lck".format(i)
        add_new_monomer(product)
        Rule('{0}_cat'.format(product), eval('{0}()'.format(previous_product)) | eval('{0}()'.format(product)),
             k_p_on_R_pmhc, k_p_off_R_pmhc)
        add_observable(product)

        no_ligand = self.unbind_ligand(i, product)

        if i == "Ls":
            no_lck = self.lck_off(no_ligand)

            Rule('{0}_uncat'.format(no_lck), eval('{0}()'.format(no_lck)) | R(), k_p_off_R, k_p_on_R)
        return product

    def cycle_3(self, i):
        if i == "Ls":
            if self.p_flag:
                Parameter('k_zap_on_R_pmhc', self.parameters['k_zap_on_R_pmhc'])
            else:
                Parameter('k_zap_on_R_pmhc', self.rate_constants.k_zap_on_R_pmhc)
            Parameter('k_zap_off_R_pmhc', self.rate_constants.k_zap_off_R_pmhc)

            # Parameter('k_lck_on_zap_R', self.rate_constants.k_lck_on_zap_R)
            Parameter('k_lck_off_zap_R', self.rate_constants.k_lck_off_zap_R)

            Parameter('k_zap_on_R', self.rate_constants.k_zap_on_R)
            Parameter('k_zap_off_R', self.rate_constants.k_zap_off_R)

        previous_product = self.cycle_2(i)
        product = "RP{0}_Lck_Zap".format(i)
        add_new_monomer(product)

        Rule('{0}_bind'.format(product), eval('{0}()'.format(previous_product)) + Zap() | eval('{0}()'.format(product)),
             k_zap_on_R_pmhc, k_zap_off_R_pmhc)
        add_observable(product)

        no_ligand = self.unbind_ligand(i, product)

        if i == "Ls":
            no_lck = self.lck_off(no_ligand, k_off="k_lck_off_zap_R")
            no_zap = self.zap_off(no_lck)

        return product

    def dephosphorylate_zap(self, product):
        new = product.replace("_Zap_P", "_Zap")
        add_new_monomer(new)
        Rule('{0}_uncat'.format(product), eval('{0}()'.format(product)) | eval('{0}()'.format(new)),
             k_p_off_R, k_p_on_R)

        return new

    def rp_zap_off(self, product):
        new = product.replace("RP_Zap_", "")
        add_new_monomer(new)
        Rule('{0}_rpzap_off'.format(product), eval('{0}()'.format(product)) | RP_Zap() + eval('{0}()'.format(new)),
             k_zap_off_R, k_zap_on_R)

        return new

    def cycle_4(self, i):

        if i == "Ls":
            if self.p_flag:
                Parameter('k_p_on_zap_species', self.parameters['k_p_on_zap_species'])
            else:
                Parameter('k_p_on_zap_species', self.rate_constants.k_p_on_zap_species)

            Parameter('k_p_off_zap_species', self.rate_constants.k_p_off_zap_species)

        previous_product = self.cycle_3(i)
        product = "RP{0}_Lck_Zap_P".format(i)
        add_new_monomer(product)

        Rule('{0}_cat'.format(product), eval('{0}()'.format(previous_product)) | eval('{0}()'.format(product)),
             k_p_on_zap_species, k_p_off_zap_species)
        add_observable(product)

        no_ligand = self.unbind_ligand(i, product)

        if i == "Ls":
            no_lck = self.lck_off(no_ligand, k_off="k_lck_off_zap_R")
            self.dephosphorylate_zap(no_lck)

        return product

    def add_negative_feedback(self):

        if self.p_flag:
            Parameter('k_neg_fb', self.parameters['k_neg_fb'])
            Parameter('k_lcki', self.parameters['k_lcki'])
        else:
            Parameter('k_neg_fb', self.rate_constants.k_negative_loop)
            Parameter('k_lcki', self.rate_constants.k_lcki)

        add_new_monomer("LckI")

        if self.self_foreign:
            Rule('Neg_fb_1_self', eval("RPLs_Lck_Zap_P()") + eval('Lck()') >> eval(
                'LckI()') + eval('RPLs_Lck_Zap_P()'), k_neg_fb)
            Rule('Neg_fb_1_foreign', eval("RPLf_Lck_Zap_P()") + eval('Lck()') >> eval(
                'LckI()') + eval("RPLf_Lck_Zap_P()"), k_neg_fb)
        else:
            Rule('Neg_fb_1', eval("RPLs_Lck_Zap_P()") + eval('Lck()') >> eval(
                'LckI()') + eval('RPLs_Lck_Zap_P()'), k_neg_fb)

        Rule('Neg_fb_2', LckI() >> Lck(), k_lcki)

    def add_negative_feedback_late_zap(self):

        if self.p_flag:
            Parameter('k_neg_fb', self.parameters['k_neg_fb'])
            Parameter('k_lcki', self.parameters['k_lcki'])
        else:
            Parameter('k_neg_fb', self.rate_constants.k_negative_loop)
            Parameter('k_lcki', self.rate_constants.k_lcki)

        add_new_monomer("ZapI")

        if self.self_foreign:
            Rule('Neg_fb_1_self', eval("RPLs_Lck_Zap_P_LATP()") + eval('Zap()') >> eval(
                'ZapI()') + eval('RPLs_Lck_Zap_P_LATP()'), k_neg_fb)
            Rule('Neg_fb_1_foreign', eval("RPLf_Lck_Zap_P_LATP()") + eval('Zap()') >> eval(
                'ZapI()') + eval("RPLf_Lck_Zap_P_LATP()"), k_neg_fb)
        else:
            Rule('Neg_fb_1', eval("RPLs_Lck_Zap_P_LATP()") + eval('Zap()') >> eval(
                'ZapI()') + eval('RPLs_Lck_Zap_P_LATP()'), k_neg_fb)

        Rule('Neg_fb_2', ZapI() >> Zap(), k_lcki)

    def add_negative_feedback_late_lat(self):

        if self.p_flag:
            Parameter('k_neg_fb', self.parameters['k_neg_fb'])
            Parameter('k_lcki', self.parameters['k_lcki'])
        else:
            Parameter('k_neg_fb', self.rate_constants.k_negative_loop)
            Parameter('k_lcki', self.rate_constants.k_lcki)

        add_new_monomer("LATI")

        if self.self_foreign:
            Rule('Neg_fb_1_self', eval("RPLs_Lck_Zap_P_LATP()") + eval('LAT()') >> eval(
                'LATI()') + eval('RPLs_Lck_Zap_P_LATP()'), k_neg_fb)
            Rule('Neg_fb_1_foreign', eval("RPLf_Lck_Zap_P_LATP()") + eval('LAT()') >> eval(
                'LATI()') + eval("RPLf_Lck_Zap_P_LATP()"), k_neg_fb)
        else:
            Rule('Neg_fb_1', eval("RPLs_Lck_Zap_P_LATP()") + eval('LAT()') >> eval(
                'LATI()') + eval('RPLs_Lck_Zap_P_LATP()'), k_neg_fb)

        Rule('Neg_fb_2', LATI() >> LAT(), k_lcki)

    def cycle_6(self, i):

        if i == "Ls":
            if self.p_flag:
                Parameter('k_lat_on_species', self.parameters['k_lat_on_species'])
            else:
                Parameter('k_lat_on_species', self.rate_constants.k_lat_on_species)
            Parameter('k_lat_off_species', self.rate_constants.k_lat_off_species)

        previous_product = self.cycle_4(i)  # self.add_negative_feedback(i)
        product = "RP{0}_Lck_Zap_P_LAT".format(i)
        add_new_monomer(product)

        Rule("{0}_bind".format(product), eval('{0}()'.format(previous_product)) + LAT() | eval('{0}()'.format(product)),
             k_lat_on_species, k_lat_off_species)
        add_observable(product)

        no_ligand = self.unbind_ligand(i, product)

        if i == "Ls":
            no_lck = self.lck_off(no_ligand, k_off="k_lck_off_zap_R")
            no_zap_p = self.dephosphorylate_zap(no_lck)
            self.rp_zap_off(no_zap_p)

        return product

    def dephosphorylate_latp(self, product):
        new = product.replace("LATP", "LAT")
        add_new_monomer(new)
        Rule('{0}_uncat'.format(product), eval('{0}()'.format(product)) | eval('{0}()'.format(new)),
             k_p_off_R, k_p_on_R)
        return new

    def cycle_7(self, i):

        if i == "Ls":
            if self.p_flag:
                Parameter('kp_on_lat1', self.parameters['kp_on_lat1'])
            else:
                Parameter('kp_on_lat1', self.rate_constants.k_p_lat_1)
            Parameter('kp_off_lat1', self.rate_constants.k_p_lat_off_species)

        previous_product = self.cycle_6(i)
        product = "RP{0}_Lck_Zap_P_LATP".format(i)
        add_new_monomer(product)

        Rule("{0}_bind".format(product), eval('{0}()'.format(previous_product)) | eval('{0}()'.format(product)),
             kp_on_lat1, kp_off_lat1)
        add_observable(product)

        no_ligand = self.unbind_ligand(i, product)

        if i == "Ls":
            no_lck = self.lck_off(no_ligand, k_off="k_lck_off_zap_R")
            no_zap_p = self.dephosphorylate_zap(no_lck)
            no_latp = self.dephosphorylate_latp(no_zap_p)

        return product

    def step_8(self, i):
        if i == "Ls":
            if self.p_flag:
                Parameter('k_p_lat_on_species', self.parameters['k_p_lat_on_species'])
            else:
                Parameter('k_p_lat_on_species', self.rate_constants.k_p_lat_2)

        previous_product = self.cycle_7(i)
        product = "LATPP"
        if i == "Ls":
            add_new_monomer(product)

        Rule("{0}_cat".format(previous_product),
             eval('{0}()'.format(previous_product)) >> eval('{0}()'.format("RP{0}_Lck_Zap_P".format(i))) + eval(
                 '{0}()'.format(product)), k_p_lat_on_species)
        if i == "Ls":
            add_observable(product)

            Rule("{0}_uncat".format(product), eval('{0}()'.format(product)) >> LAT(), k_p_off_R_pmhc)

        return product

    def non_specific_step_8(self, i):

        if i == "Ls":
            Parameter('k_product_on', self.rate_constants.k_product_on)

        previous_product = self.step_8(i)
        product = "final_product"
        if i == "Ls":
            add_new_monomer(product)
            Rule("{0}_convert".format(product), eval('{0}()'.format(previous_product)) | eval('{0}()'.format(product)),
                 k_product_on, k_p_off_R_pmhc)
            add_observable(product)

        return product

    def add_step_8_sos(self, i):

        previous_product = self.step_8(i)
        product = previous_product + "_Sos"

        if i == "Ls":
            if self.p_flag:
                Parameter('k_sos_on', self.parameters['k_8_1'])
            else:
                Parameter('k_sos_on', self.rate_constants.k_sos_on)

            Parameter('k_sos_off', self.rate_constants.k_sos_off)

            add_new_monomer(product)
            Rule("{0}_bind".format(product),
                 eval('{0}()'.format(previous_product)) + Sos() | eval('{0}()'.format(product)),
                 k_sos_on, k_sos_off)
            add_observable(product)

        return product

    def add_step_9(self, i):
        previous_product = self.add_step_8_sos(i)
        product_rd = previous_product + "_Ras_GDP"

        if i == "Ls":
            Parameter('k_sos_on_rgdp', self.rate_constants.k_sos_on_rgdp)
            Parameter('k_sos_off_rgdp', self.rate_constants.k_sos_off_rgdp)

            add_new_monomer(product_rd)

            Rule('{0}_bind'.format(product_rd),
                 eval('{0}()'.format(previous_product)) + Ras_GDP() | eval('{0}()'.format(product_rd)),
                 k_sos_on_rgdp, k_sos_off_rgdp)

        product_rt = previous_product + "_Ras_GTP"

        if i == "Ls":
            Parameter('k_sos_on_rgtp', self.rate_constants.k_sos_on_rgtp)
            Parameter('k_sos_off_rgtp', self.rate_constants.k_sos_off_rgtp)

            add_new_monomer(product_rt)

            Rule('{0}_bind'.format(product_rt),
                 eval('{0}()'.format(previous_product)) + Ras_GTP() | eval('{0}()'.format(product_rt)),
                 k_sos_on_rgtp, k_sos_off_rgtp)

        return product_rd, product_rt

    def add_step_10(self, i):
        previous_product_rd, previous_product_rt = self.add_step_9(i)
        product_rt_rd = previous_product_rt + "_Ras_GDP"

        if i == "Ls":
            Parameter('k_rgdp_on_sos_rgtp', self.rate_constants.k_rgdp_on_sos_rgtp)
            Parameter('k_rgdp_off_sos_rgtp', self.rate_constants.k_rgdp_off_sos_rgtp)
            Parameter('k_cat_3', self.rate_constants.k_cat_3)

            add_new_monomer(product_rt_rd)

            Rule('{0}_bind'.format(product_rt_rd),
                 eval('{0}()'.format(previous_product_rt)) + Ras_GDP() | eval('{0}()'.format(product_rt_rd)),
                 k_rgdp_on_sos_rgtp, k_rgdp_off_sos_rgtp)

            Rule('{0}_cat'.format(product_rt_rd),
                 eval('{0}()'.format(product_rt_rd)) >> eval('{0}()'.format(previous_product_rt)) + Ras_GTP(),
                 k_cat_3)

        product_rd_rd = previous_product_rd + "_Ras_GDP"

        if i == "Ls":
            Parameter('k_rgdp_on_sos_rgdp', self.rate_constants.k_rgdp_on_sos_rgdp)
            Parameter('k_rgdp_off_sos_rgdp', self.rate_constants.k_rgdp_off_sos_rgdp)
            Parameter('k_cat_4', self.rate_constants.k_cat_4)

            add_new_monomer(product_rd_rd)

            Rule('{0}_bind'.format(product_rd_rd),
                 eval('{0}()'.format(previous_product_rd)) + Ras_GDP() | eval('{0}()'.format(product_rd_rd)),
                 k_rgdp_on_sos_rgdp, k_rgdp_off_sos_rgdp)

            Rule('{0}_cat'.format(product_rd_rd),
                 eval('{0}()'.format(product_rd_rd)) >> eval('{0}()'.format(previous_product_rd)) + Ras_GTP(),
                 k_cat_4)

        # Deactivate - convert Ras_GTP to Ras_GDP

        product = "Ras_GAP_Ras_GTP"
        new_product = "Ras_GTP"

        if i == "Ls":
            Parameter('k_rgap_on_rgtp', self.rate_constants.k_rgap_on_rgtp)
            Parameter('k_rgap_off_rgtp', self.rate_constants.k_rgap_off_rgtp)
            Parameter('k_cat_5', self.rate_constants.k_cat_5)

            add_new_monomer(product)

            Rule('{0}_bind'.format(product),
                 Ras_GAP() + Ras_GTP() | eval('{0}()'.format(product)),
                 k_rgap_on_rgtp, k_rgap_off_rgtp)

            add_observable(product)

            Rule('{0}_cat'.format(product),
                 eval('{0}()'.format(product)) >> Ras_GAP() + Ras_GDP(),
                 k_cat_5)

            add_observable(new_product)

        return new_product

    def add_positive_feedback_nonspecific(self, i):

        if i == "Ls":
            Parameter('k_positive_fb', self.rate_constants.k_positive_loop)

        product = self.non_specific_step_8(i)
        if i == "Ls":
            Rule("{0}_posfb".format(product),
                 eval('{0}()'.format(product)) + eval('{0}()'.format("LATP")) >> eval(
                     '{0}()'.format(product)) + eval('{0}()'.format(product)), k_positive_fb)

        return product

    def make_model(self):

        self.define_monomers()
        self.define_rates()

        observables = []
        for i in self.output:
            if self.steps == 0:
                product = self.add_step_0(i)
            elif self.steps == 1:
                product = self.cycle_1(i)
            elif self.steps == 2:
                product = self.cycle_2(i)
            elif self.steps == 3:
                product = self.cycle_3(i)
            elif self.steps == 4:
                product = self.cycle_4(i)
            elif self.steps == 5:
                product = self.cycle_4(i)
            elif self.steps == 6:
                product = self.cycle_6(i)
            elif self.steps == 7:
                product = self.cycle_7(i)
            elif self.steps == 8:
                product = self.step_8(i)
            elif self.steps == 9:
                product = self.add_step_8_sos(i)
            else:
                product = self.add_step_10(i)

            if "O_{0}".format(product) not in observables:
                observables.append("O_{0}".format(product))

        # NEGATIVE FEEDBACK LOOP has been added here!!!!

        if self.steps > 4:
            self.add_negative_feedback()

        # if self.steps > 7:
        #     self.add_negative_feedback_late_lat()

        return observables

    def main_truncated_time(self):
        output = []
        ligand_array = []
        observables = self.make_model()

        write_columns(observables)
        write_model_attributes(self.model.rules, "rules")
        write_model_attributes(self.model.parameters, "parameters")
        write_model_attributes(self.model.observables, "observables")

        np.savetxt("time", self.tspan, fmt='%f')
        time_index = 2

        for ligand in self.p_ligand:
            self.model.parameters['Ls_0'].value = ligand
            ligand_array.append(ligand)

            y = odesolve(self.model, self.tspan, compiler="python")

            if len(observables) > 1:
                output_array = y[observables[0]] + y[observables[1]]
                if self.num_samples == 1:
                    np.savetxt("{0}_output".format(observables[0]), y[observables[0]], fmt='%f')
                    np.savetxt("{0}_output".format(observables[1]), y[observables[1]], fmt='%f')
                    np.savetxt("output_array", output_array, fmt='%f')

                output.append(output_array[time_index])
            else:
                output_array = y[observables[0]]
                output.append(output_array[time_index])

        np.savetxt("truncated_time", [self.tspan[time_index]], fmt='%f')
        np.savetxt("Ligand_concentrations", ligand_array, fmt='%f')
        np.savetxt("output", output, fmt='%f')

    def main(self):
        output = []
        ligand_array = []
        observables = self.make_model()

        write_columns(observables)
        write_model_attributes(self.model.rules, "rules")
        write_model_attributes(self.model.parameters, "parameters")
        write_model_attributes(self.model.observables, "observables")

        np.savetxt("time", self.tspan, fmt='%f')

        for ligand in self.p_ligand:
            self.model.parameters['Ls_0'].value = ligand
            ligand_array.append(ligand)

            y = odesolve(self.model, self.tspan, compiler="python")

            if len(observables) > 1:
                output_array = y[observables[0]] + y[observables[1]]
                if self.num_samples == 1:
                    np.savetxt("{0}_output".format(observables[0]), y[observables[0]], fmt='%f')
                    np.savetxt("{0}_output".format(observables[1]), y[observables[1]], fmt='%f')
                    np.savetxt("output_array", output_array, fmt='%f')

                output.append(output_array[-1])
            else:
                output_array = y[observables[0]]
                output.append(output_array[-1])

        np.savetxt("Ligand_concentrations", ligand_array, fmt='%f')
        np.savetxt("output", output, fmt='%f')


class NonSpecificEarlyPositiveFeedback(PysbTcrSelfWithForeign):
    def __init__(self, steps=3, self_foreign=False, lf=30):
        PysbTcrSelfWithForeign.__init__(self, steps=steps, self_foreign=self_foreign, lf=lf)

    def step_4(self, i):

        if i == "Ls":
            if self.p_flag:
                Parameter('k_p_on_zap_species', self.parameters['k_p_on_zap_species'])
            else:
                Parameter('k_p_on_zap_species', self.rate_constants.k_p_on_zap_species)

            Parameter('k_p_off_zap_species', self.rate_constants.k_p_off_zap_species)

        previous_product = self.cycle_3(i)
        product = "Zap_P"
        if i == "Ls":
            add_new_monomer(product)

        Rule('{0}_cat'.format(previous_product),
             eval('{0}()'.format(previous_product)) >> eval('{0}()'.format("RP{0}_Lck".format(i))) +
             eval('{0}()'.format(product)), k_p_on_zap_species)

        if i == "Ls":
            add_observable(product)

            Rule("{0}_uncat".format(product), eval('{0}()'.format(product)) >> Zap(), k_p_off_R_pmhc)

        return product

    def add_step_9(self, i):
        previous_product = self.step_4(i)
        product_rd = previous_product + "_Ras_GDP"

        if i == "Ls":
            Parameter('k_sos_on_rgdp', self.rate_constants.k_sos_on_rgdp)
            Parameter('k_sos_off_rgdp', self.rate_constants.k_sos_off_rgdp)

            add_new_monomer(product_rd)

            Rule('{0}_bind'.format(product_rd),
                 eval('{0}()'.format(previous_product)) + Ras_GDP() | eval('{0}()'.format(product_rd)),
                 k_sos_on_rgdp, k_sos_off_rgdp)

        product_rt = previous_product + "_Ras_GTP"

        if i == "Ls":
            Parameter('k_sos_on_rgtp', self.rate_constants.k_sos_on_rgtp)
            Parameter('k_sos_off_rgtp', self.rate_constants.k_sos_off_rgtp)

            add_new_monomer(product_rt)

            Rule('{0}_bind'.format(product_rt),
                 eval('{0}()'.format(previous_product)) + Ras_GTP() | eval('{0}()'.format(product_rt)),
                 k_sos_on_rgtp, k_sos_off_rgtp)

        return product_rd, product_rt

    def make_model(self):

        self.define_monomers()
        self.define_rates()

        observables = []
        for i in self.output:
            if self.steps == 0:
                product = self.add_step_0(i)
            elif self.steps == 1:
                product = self.cycle_1(i)
            elif self.steps == 2:
                product = self.cycle_2(i)
            elif self.steps == 3:
                product = self.cycle_3(i)
            elif self.steps == 4:
                product = self.step_4(i)
            else:
                product = self.add_step_10(i)

            if "O_{0}".format(product) not in observables:
                observables.append("O_{0}".format(product))

        return observables


class EarlyPositiveFeedback(PysbTcrSelfWithForeign):
    def __init__(self, steps=3, self_foreign=False, lf=30):
        PysbTcrSelfWithForeign.__init__(self, steps=steps, self_foreign=self_foreign, lf=lf)

    def add_step_8_sos(self, i):
        if i == "Ls":
            if self.p_flag:
                Parameter('k_sos_on', self.parameters['k_8_1'])
            else:
                Parameter('k_sos_on', self.rate_constants.k_sos_on)

            Parameter('k_sos_off', self.rate_constants.k_sos_off)

        previous_product = self.cycle_3(i)
        product = previous_product + "_Sos"

        add_new_monomer(product)
        Rule("{0}_bind".format(product), eval('{0}()'.format(previous_product)) + Sos() | eval('{0}()'.format(product)),
             k_sos_on, k_sos_off)
        add_observable(product)

        return product

    def add_step_9(self, i):
        if i == "Ls":
            Parameter('k_sos_on_rgdp', self.rate_constants.k_sos_on_rgdp)
            Parameter('k_sos_off_rgdp', self.rate_constants.k_sos_off_rgdp)

        previous_product = self.cycle_3(i)
        product_rd = previous_product + "_Ras_GDP"

        add_new_monomer(product_rd)

        Rule('{0}_bind'.format(product_rd),
             eval('{0}()'.format(previous_product)) + Ras_GDP() | eval('{0}()'.format(product_rd)),
             k_sos_on_rgdp, k_sos_off_rgdp)

        product_rt = previous_product + "_Ras_GTP"

        if i == "Ls":
            Parameter('k_sos_on_rgtp', self.rate_constants.k_sos_on_rgtp)
            Parameter('k_sos_off_rgtp', self.rate_constants.k_sos_off_rgtp)

        add_new_monomer(product_rt)

        Rule('{0}_bind'.format(product_rt),
             eval('{0}()'.format(previous_product)) + Ras_GTP() | eval('{0}()'.format(product_rt)),
             k_sos_on_rgtp, k_sos_off_rgtp)

        return product_rd, product_rt

    def add_step_10(self, i):
        previous_product_rd, previous_product_rt = self.add_step_9(i)
        product_rt_rd = previous_product_rt + "_Ras_GDP"

        if i == "Ls":
            Parameter('k_rgdp_on_sos_rgtp', self.rate_constants.k_rgdp_on_sos_rgtp)
            Parameter('k_rgdp_off_sos_rgtp', self.rate_constants.k_rgdp_off_sos_rgtp)
            Parameter('k_cat_3', self.rate_constants.k_cat_3)

        add_new_monomer(product_rt_rd)

        Rule('{0}_bind'.format(product_rt_rd),
             eval('{0}()'.format(previous_product_rt)) + Ras_GDP() | eval('{0}()'.format(product_rt_rd)),
             k_rgdp_on_sos_rgtp, k_rgdp_off_sos_rgtp)

        Rule('{0}_cat'.format(product_rt_rd),
             eval('{0}()'.format(product_rt_rd)) >> eval('{0}()'.format(previous_product_rt)) + Ras_GTP(),
             k_cat_3)

        product_rd_rd = previous_product_rd + "_Ras_GDP"

        if i == "Ls":
            Parameter('k_rgdp_on_sos_rgdp', self.rate_constants.k_rgdp_on_sos_rgdp)
            Parameter('k_rgdp_off_sos_rgdp', self.rate_constants.k_rgdp_off_sos_rgdp)
            Parameter('k_cat_4', self.rate_constants.k_cat_4)

        add_new_monomer(product_rd_rd)

        Rule('{0}_bind'.format(product_rd_rd),
             eval('{0}()'.format(previous_product_rd)) + Ras_GDP() | eval('{0}()'.format(product_rd_rd)),
             k_rgdp_on_sos_rgdp, k_rgdp_off_sos_rgdp)

        Rule('{0}_cat'.format(product_rd_rd),
             eval('{0}()'.format(product_rd_rd)) >> eval('{0}()'.format(previous_product_rd)) + Ras_GTP(),
             k_cat_4)

        # Deactivate - convert Ras_GTP to Ras_GDP

        product = "Ras_GAP_Ras_GTP"
        new_product = "Ras_GTP"

        if i == "Ls":
            Parameter('k_rgap_on_rgtp', self.rate_constants.k_rgap_on_rgtp)
            Parameter('k_rgap_off_rgtp', self.rate_constants.k_rgap_off_rgtp)
            Parameter('k_cat_5', self.rate_constants.k_cat_5)

            add_new_monomer(product)

            Rule('{0}_bind'.format(product),
                 Ras_GAP() + Ras_GTP() | eval('{0}()'.format(product)),
                 k_rgap_on_rgtp, k_rgap_off_rgtp)

            add_observable(product)

            Rule('{0}_cat'.format(product),
                 eval('{0}()'.format(product)) >> Ras_GAP() + Ras_GDP(),
                 k_cat_5)

            add_observable(new_product)

        return new_product

    def make_model(self):

        self.define_monomers()
        self.define_rates()

        observables = []
        for i in self.output:
            if self.steps == 0:
                product = self.add_step_0(i)
            elif self.steps == 1:
                product = self.cycle_1(i)
            elif self.steps == 2:
                product = self.cycle_2(i)
            elif self.steps == 3:
                product = self.cycle_3(i)
            else:
                product = self.add_step_10(i)

            if "O_{0}".format(product) not in observables:
                observables.append("O_{0}".format(product))

        return observables


class LatPhosphorylationExtension(PysbTcrSelfWithForeign):
    def __init__(self, steps=8, self_foreign=False, lf=30):
        PysbTcrSelfWithForeign.__init__(self, steps=steps, self_foreign=self_foreign, lf=lf)

    def cycle_7(self, i):

        if i == "Ls":
            if self.p_flag:
                Parameter('kp_on_lat1', self.parameters['kp_on_lat1'])
            else:
                Parameter('kp_on_lat1', self.rate_constants.k_p_lat_2)
            Parameter('kp_off_lat1', self.rate_constants.k_p_lat_off_species)

        previous_product = self.cycle_6(i)
        product = "RP{0}_Lck_Zap_P_LATP".format(i)
        add_new_monomer(product)

        Rule("{0}_bind".format(product), eval('{0}()'.format(previous_product)) | eval('{0}()'.format(product)),
             kp_on_lat1, kp_off_lat1)
        add_observable(product)

        no_ligand = self.unbind_ligand(i, product)

        if i == "Ls":
            no_lck = self.lck_off(no_ligand, k_off="k_lck_off_zap_R")
            no_zap_p = self.dephosphorylate_zap(no_lck)
            no_latp = self.dephosphorylate_latp(no_zap_p)

        return product

    def cycle_8(self, i):

        if i == "Ls":
            if self.p_flag:
                Parameter('kp_on_lat2', self.parameters['kp_on_lat2'])
            else:
                Parameter('kp_on_lat2', self.rate_constants.k_p_lat_2)
            Parameter('kp_off_lat2', self.rate_constants.k_p_lat_off_species)

        previous_product = self.cycle_7(i)
        product = "RP{0}_Lck_Zap_P_LATPP".format(i)
        add_new_monomer(product)

        Rule("{0}_bind".format(product), eval('{0}()'.format(previous_product)) | eval('{0}()'.format(product)),
             kp_on_lat2, kp_off_lat2)
        add_observable(product)

        no_ligand = self.unbind_ligand(i, product)

        if i == "Ls":
            no_lck = self.lck_off(no_ligand, k_off="k_lck_off_zap_R")
            no_zap_p = self.dephosphorylate_zap(no_lck)
            no_latp = self.dephosphorylate_latp(no_zap_p)

        return product

    def make_model(self):

        self.define_monomers()
        self.define_rates()

        observables = []
        for i in self.output:
            if self.steps == 0:
                product = self.add_step_0(i)
            elif self.steps == 1:
                product = self.cycle_1(i)
            elif self.steps == 2:
                product = self.cycle_2(i)
            elif self.steps == 3:
                product = self.cycle_3(i)
            elif self.steps == 4:
                product = self.cycle_4(i)
            elif self.steps == 5:
                product = self.cycle_4(i)
            elif self.steps == 6:
                product = self.cycle_6(i)
            elif self.steps == 7:
                product = self.cycle_7(i)
            else:
                product = self.cycle_8(i)

            if "O_{0}".format(product) not in observables:
                observables.append("O_{0}".format(product))

        return observables


class PysbTcrLckFeedbackLoop(PysbTcrSelfWithForeign):
    def __init__(self, steps=8, self_foreign=False, lf=30):
        PysbTcrSelfWithForeign.__init__(self, steps=steps, self_foreign=self_foreign, lf=lf)

    def inactive_lck_off(self, no_ligand_product, k_off="k_lck_off_R"):
        new = no_ligand_product.replace("_LckI", "")
        self.add_new_monomer(new)
        Rule('{0}_unbindLckI'.format(no_ligand_product),
             eval('{0}()'.format(no_ligand_product)) | eval('{0}()'.format(new)) + LckI(),
             eval(k_off), k_lck_on_R)
        return new

    def add_negative_feedback(self):

        if self.p_flag:
            Parameter('k_lcki', self.parameters['k_lcki'])
        else:
            Parameter('k_lcki', self.rate_constants.k_lcki)

        Rule('Neg_fb_2', LckI() >> Lck(), k_lcki)

    def cycle_4(self, i):

        if i == "Ls":
            if self.p_flag:
                Parameter('k_p_on_zap_species', self.parameters['k_p_on_zap_species'])
            else:
                Parameter('k_p_on_zap_species', self.rate_constants.k_p_on_zap_species)

            Parameter('k_p_off_zap_species', self.rate_constants.k_p_off_zap_species)

        previous_product = self.cycle_3(i)
        product = "RP{0}_LckI_Zap_P".format(i)
        add_new_monomer(product)

        Rule('{0}_cat'.format(product), eval('{0}()'.format(previous_product)) | eval('{0}()'.format(product)),
             k_p_on_zap_species, k_p_off_zap_species)
        add_observable(product)

        no_ligand = self.unbind_ligand(i, product)

        if i == "Ls":
            no_lck = self.inactive_lck_off(no_ligand, k_off="k_lck_off_zap_R")
            self.dephosphorylate_zap(no_lck)

        return product

    def cycle_6(self, i):

        if i == "Ls":
            if self.p_flag:
                Parameter('k_lat_on_species', self.parameters['k_lat_on_species'])
            else:
                Parameter('k_lat_on_species', self.rate_constants.k_lat_on_species)
            Parameter('k_lat_off_species', self.rate_constants.k_lat_off_species)

            add_new_monomer("LckI")

        previous_product = self.cycle_4(i)
        product = "RP{0}_LckI_Zap_P_LAT".format(i)
        add_new_monomer(product)

        Rule("{0}_bind".format(product), eval('{0}()'.format(previous_product)) + LAT() | eval('{0}()'.format(product)),
             k_lat_on_species, k_lat_off_species)
        add_observable(product)

        no_ligand = self.unbind_ligand(i, product)

        if i == "Ls":
            no_lck = self.inactive_lck_off(no_ligand, k_off="k_lck_off_zap_R")
            no_zap_p = self.dephosphorylate_zap(no_lck)
            self.rp_zap_off(no_zap_p)

        return product

    def cycle_7(self, i):

        if i == "Ls":
            if self.p_flag:
                Parameter('kp_on_lat1', self.parameters['kp_on_lat1'])
            else:
                Parameter('kp_on_lat1', self.rate_constants.k_p_lat_1)
            Parameter('kp_off_lat1', self.rate_constants.k_p_lat_off_species)

        previous_product = self.cycle_6(i)
        product = "RP{0}_LckI_Zap_P_LATP".format(i)
        add_new_monomer(product)

        Rule("{0}_bind".format(product), eval('{0}()'.format(previous_product)) | eval('{0}()'.format(product)),
             kp_on_lat1, kp_off_lat1)
        add_observable(product)

        no_ligand = self.unbind_ligand(i, product)

        if i == "Ls":
            no_lck = self.inactive_lck_off(no_ligand, k_off="k_lck_off_zap_R")
            no_zap_p = self.dephosphorylate_zap(no_lck)
            no_latp = self.dephosphorylate_latp(no_zap_p)

        return product

    def step_8(self, i):
        if i == "Ls":
            if self.p_flag:
                Parameter('k_p_lat_on_species', self.parameters['k_p_lat_on_species'])
            else:
                Parameter('k_p_lat_on_species', self.rate_constants.k_p_lat_2)

        previous_product = self.cycle_7(i)
        product = "LATPP"
        if i == "Ls":
            add_new_monomer(product)

        Rule("{0}_cat".format(previous_product),
             eval('{0}()'.format(previous_product)) >> eval('{0}()'.format("RP{0}_LckI_Zap_P".format(i))) + eval(
                 '{0}()'.format(product)), k_p_lat_on_species)
        if i == "Ls":
            add_observable(product)

            Rule("{0}_uncat".format(product), eval('{0}()'.format(product)) >> LAT(), k_p_off_R_pmhc)

        return product

    def make_model(self):

        self.define_monomers()
        self.define_rates()

        observables = []
        for i in self.output:
            if self.steps == 0:
                product = self.add_step_0(i)
            elif self.steps == 1:
                product = self.cycle_1(i)
            elif self.steps == 2:
                product = self.cycle_2(i)
            elif self.steps == 3:
                product = self.cycle_3(i)
            elif self.steps == 4:
                product = self.cycle_4(i)
            elif self.steps == 5:
                product = self.cycle_4(i)
            elif self.steps == 6:
                product = self.cycle_6(i)
            elif self.steps == 7:
                product = self.cycle_7(i)
            else:
                product = self.step_8(i)

            if "O_{0}".format(product) not in observables:
                observables.append("O_{0}".format(product))

        # NEGATIVE FEEDBACK LOOP has been added here!!!!

        if self.steps > 3:
            self.add_negative_feedback()

        return observables


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Submitting ode calculations as function of steps",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--steps', dest='steps', action='store', type=int, help="number of KP steps.")
    parser.add_argument('--lf', dest='lf', action='store', type=int, help="Number of foreign ligands.")
    parser.add_argument('--early_pos_fb', dest='early_pos_fb', action='store_true', default=False,
                        help='Flag for building and submitting early positive feedback loop.')

    args = parser.parse_args()

    if args.lf:
        if args.early_pos_fb:
            tcr = EarlyPositiveFeedback(steps=args.steps, self_foreign=True, lf=args.lf)

        else:
            tcr = PysbTcrSelfWithForeign(steps=args.steps, self_foreign=True, lf=args.lf)

    else:
        if args.early_pos_fb:
            tcr = EarlyPositiveFeedback(steps=args.steps)

        else:
            tcr = PysbTcrSelfWithForeign(steps=args.steps)

    tcr.main()

## Uncomment to make reaction network

# # tcr = PysbTcrSelfWithForeign(steps=8, self_foreign=True)
# tcr = PysbTcrLckFeedbackLoop(steps=8, self_foreign=True)
#
# # tcr = LatPhosphorylationExtension(steps=7)
# tcr.make_model()

# tcr = PysbTcrCubicFbLoop(steps=8, self_foreign=True)
# tcr.make_model()
# tcr = PysbTcrCubicFbLoop(steps=8, self_foreign=True)
# tcr.make_model()

## python /anaconda2/lib/python2.7/site-packages/pysb/tools/render_reactions.py ../SSC_python_modules/pysb_t_cell_network.py > model.dot
## dot -Tpdf model.dot -o model.pdf
