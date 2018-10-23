import os
import pickle

import numpy as np
from pysb import *
from pysb.integrate import odesolve
from pysb.simulator.bng import BngSimulator

from realistic_network import BindingParameters, InitialConcentrations


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

        self.lf = lf
        self.ls = 1000
        self.steps = steps

        self.run_time = 10000
        self.tspan = np.linspace(0, self.run_time)

        self.mu = 6
        self.sigma = 1.0
        self.num_samples = 1000

        self.p_flag = os.path.exists("parameters.pickle")

        if self.p_flag:
            self.parameters = pickle.load(open("parameters.pickle", "rb"))
            print(self.parameters)

        self.p_ligand = [int(i) for i in np.round(np.random.lognormal(self.mu, self.sigma, self.num_samples))]

        if self_foreign:
            self.output = ["Ls", "Lf"]
        else:
            self.output = ["Ls"]

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

    def define_rates(self):
        Parameter('k_L_on', self.rate_constants.k_L_on)
        if "Lf" in self.output:
            Parameter('k_Lf_off', self.rate_constants.k_foreign_off)
        Parameter('k_Ls_off', self.rate_constants.k_self_off)

    def add_step_0(self, i):
        product = "R{0}".format(i)
        add_new_monomer(product)

        Rule('{0}_bind'.format(product), R() + eval('{0}()'.format(i)) | eval('{0}()'.format(product)),
             k_L_on, eval('k_{0}_off'.format(i)))

        add_observable(i)
        add_observable(product)
        return product

    def zap_off(self, product):
        new = product.replace("_Zap", "")
        add_new_monomer(new)
        Rule('{0}_unbindZap'.format(product), eval('{0}()'.format(product)) | eval('{0}()'.format(new)) + Zap(),
             k_zap_off_R, k_zap_on_R)
        return new

    def lck_off(self, no_ligand_product, k_off="k_lck_off_R"):
        new = no_ligand_product.replace("_Lck", "")
        add_new_monomer(new)
        Rule('{0}_unbindLck'.format(no_ligand_product),
             eval('{0}()'.format(no_ligand_product)) | eval('{0}()'.format(new)) + Lck(),
             eval(k_off), k_lck_on_R)

        return new

    def unbind_ligand(self, i, final_product):
        new = final_product.replace(i, "")
        add_new_monomer(new)
        Rule('{0}_unbindL'.format(final_product),
             eval('{0}()'.format(final_product)) | eval('{0}()'.format(i)) + eval('{0}()'.format(new)),
             eval('k_{0}_off'.format(i)), k_L_on)

        return new

    def cycle_1(self, i):
        if i == "Ls":
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

    def add_negative_feedback(self, i):

        if i == "Ls":
            Parameter('k_negative_fb', self.rate_constants.k_negative_loop)
            Parameter('k_diss_rllcki', 20.0)
            Parameter('k_lcki', 20.0)

        previous_product = self.cycle_4(i)
        add_new_monomer("R{0}_LckI".format(i))

        Rule('{0}_negfb_1'.format(previous_product),
             eval('{0}()'.format(previous_product)) + eval('R{0}_Lck()'.format(i)) >> eval(
                 'R{0}_LckI()'.format(i)) + eval('{0}()'.format(previous_product)),
             k_negative_fb)

        add_new_monomer("LckI")

        Rule('{0}_negfb_2'.format('R{0}_LckI'.format(i)),
             eval('R{0}_LckI()'.format(i)) >> eval('R{0}()'.format(i)) + eval('LckI()'),
             k_diss_rllcki)

        if i == "Ls":
            Rule('{0}_negfb_3'.format("LckI"), LckI() >> Lck(), k_lcki)

        return previous_product

    def cycle_6(self, i):

        if i == "Ls":
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

    def non_specific_step_7(self, i):
        if i == "Ls":
            # if self.p_flag:
            #     Parameter('k_p_lat_on_species', self.parameters['k_7'])
            # else:
            Parameter('k_p_lat_on_species', self.rate_constants.k_p_lat_on_species)

        previous_product = self.cycle_6(i)
        product = "LATP"
        if i == "Ls":
            add_new_monomer(product)

        Rule("{0}_cat".format(previous_product),
             eval('{0}()'.format(previous_product)) >> eval('{0}()'.format("RP{0}_Lck_Zap_P".format(i))) + eval(
                 '{0}()'.format(product)), k_p_lat_on_species)
        if i == "Ls":
            add_observable(product)

            Rule("{0}_uncat".format(product), eval('{0}()'.format(product)) >> LAT(), k_p_off_R_pmhc)

        return product

    # def add_step_7(self, i):
    #
    #     if i == "Ls":
    #         Parameter('k_p_lat_on_species', self.rate_constants.k_p_lat_on_species)
    #
    #     previous_product = self.cycle_6(i)
    #     product = "{0}_LATP".format(i)
    #     self.add_new_monomer(product)
    #
    #     Rule("{0}_cat".format(product),
    #          eval('{0}()'.format(previous_product)) >> eval('{0}()'.format("RP{0}_Lck_Zap_P".format(i))) + eval(
    #              '{0}()'.format(product)), k_p_lat_on_species)
    #     self.add_observable(product)
    #
    #     Rule("{0}_uncat".format(product), eval('{0}()'.format(product)) >> LAT(), k_p_off_R_pmhc)
    #
    #     return product

    def non_specific_step_8(self, i):

        if i == "Ls":
            Parameter('k_product_on', self.rate_constants.k_product_on)

        previous_product = self.non_specific_step_7(i)
        product = "final_product"
        if i == "Ls":
            add_new_monomer(product)
            Rule("{0}_convert".format(product), eval('{0}()'.format(previous_product)) | eval('{0}()'.format(product)),
                 k_product_on, k_p_off_R_pmhc)
            add_observable(product)

        return product

    def add_step_8_sos(self, i):

        previous_product = self.non_specific_step_7(i)
        product = "LATP_Sos"

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
            Parameter('k_sos_on_rgdp', 0.0024)
            Parameter('k_sos_off_rgdp', 3.0)

            add_new_monomer(product_rd)

            Rule('{0}_bind'.format(product_rd),
                 eval('{0}()'.format(previous_product)) + Ras_GDP() | eval('{0}()'.format(product_rd)),
                 k_sos_on_rgdp, k_sos_off_rgdp)

        product_rt = previous_product + "_Ras_GTP"

        if i == "Ls":
            Parameter('k_sos_on_rgtp', 0.0022)
            Parameter('k_sos_off_rgtp', 0.4)

            add_new_monomer(product_rt)

            Rule('{0}_bind'.format(product_rt),
                 eval('{0}()'.format(previous_product)) + Ras_GTP() | eval('{0}()'.format(product_rt)),
                 k_sos_on_rgtp, k_sos_off_rgtp)

        return product_rd, product_rt

    def add_step_10(self, i):
        previous_product_rd, previous_product_rt = self.add_step_9(i)
        product_rt_rd = previous_product_rt + "_Ras_GDP"

        if i == "Ls":
            Parameter('k_rgdp_on_sos_rgtp', 0.001)
            Parameter('k_rgdp_off_sos_rgtp', 0.1)
            Parameter('k_cat_3', 0.038 * 2.5)

            add_new_monomer(product_rt_rd)

            Rule('{0}_bind'.format(product_rt_rd),
                 eval('{0}()'.format(previous_product_rt)) + Ras_GDP() | eval('{0}()'.format(product_rt_rd)),
                 k_rgdp_on_sos_rgtp, k_rgdp_off_sos_rgtp)

            Rule('{0}_cat'.format(product_rt_rd),
                 eval('{0}()'.format(product_rt_rd)) >> eval('{0}()'.format(previous_product_rt)) + Ras_GTP(),
                 k_cat_3)

        product_rd_rd = previous_product_rd + "_Ras_GDP"

        if i == "Ls":
            Parameter('k_rgdp_on_sos_rgdp', 0.0014)
            Parameter('k_rgdp_off_sos_rgdp', 1.0)
            Parameter('k_cat_4', 0.003)

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
            Parameter('k_rgap_on_rgtp', 0.0348)
            Parameter('k_rgap_off_rgtp', 0.2)
            Parameter('k_cat_5', 0.1)

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

    # def add_step_9_sos(self, i):
    #
    #     previous_product = self.add_step_8_sos(i)
    #     product = "LATP_Grb_Sos"
    #
    #     if i == "Ls":
    #         if self.p_flag:
    #             Parameter('k_sos_on', self.parameters['k_9_1'])
    #         else:
    #             Parameter('k_sos_on', 0.0002)
    #
    #         Parameter('k_sos_off', 0.005)
    #
    #         add_new_monomer(product)
    #         Rule("{0}_bind".format(product), eval('{0}()'.format(previous_product)) + Sos() | eval('{0}()'.format(product)),
    #              k_sos_on, k_sos_off)
    #         add_observable(product)
    #
    #     return product

    # def add_step_8(self, i):
    #
    #     if i == "Ls":
    #         Parameter('k_product_on', self.rate_constants.k_product_on)
    #
    #     previous_product = self.add_step_7(i)
    #     product = "{0}_product".format(i)
    #     self.add_new_monomer(product)
    #
    #     Rule("{0}_convert".format(product), eval('{0}()'.format(previous_product)) | eval('{0}()'.format(product)),
    #          k_product_on, k_p_off_R_pmhc)
    #     self.add_observable(product)
    #
    #     return product

    def add_positive_feedback_nonspecific(self, i):

        if i == "Ls":
            Parameter('k_positive_fb', self.rate_constants.k_positive_loop)

        product = self.non_specific_step_8(i)
        if i == "Ls":
            Rule("{0}_posfb".format(product),
                 eval('{0}()'.format(product)) + eval('{0}()'.format("LATP")) >> eval(
                     '{0}()'.format(product)) + eval('{0}()'.format(product)), k_positive_fb)

        return product

    # def add_positive_feedback(self, i):
    #
    #     if i == "Ls":
    #         Parameter('k_positive_fb', self.rate_constants.k_positive_loop)
    #
    #     product = self.add_step_8(i)
    #
    #     Rule("{0}_posfb".format(product),
    #          eval('{0}()'.format(product)) + eval('{0}()'.format("{0}_LATP".format(i))) >> eval(
    #              '{0}()'.format(product)) + eval('{0}()'.format(product)), k_positive_fb)
    #
    #     return product

    def make_model(self):
        Model()

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
                product = self.add_negative_feedback(i)
            elif self.steps == 6:
                product = self.cycle_6(i)
            elif self.steps == 7:
                product = self.non_specific_step_7(i)
            elif self.steps == 8:
                # product = self.non_specific_step_8(i)
                product = self.add_step_8_sos(i)
            else:
                product = self.add_step_10(i)

                # product = self.add_step_9_sos(i)
                # product = self.non_specific_step_8(i)
            # else:
            #     # product = self.add_positive_feedback(i)
            #     product = self.add_positive_feedback_nonspecific(i)

            if "O_{0}".format(product) not in observables:
                observables.append("O_{0}".format(product))

        return observables

    def main(self):
        output = []
        ligand_array = []
        observables = self.make_model()
        # print(observables)

        write_columns(observables)
        write_model_attributes(model.rules, "rules")
        write_model_attributes(model.parameters, "parameters")
        write_model_attributes(model.observables, "observables")

        np.savetxt("time", self.tspan, fmt='%f')

        for ligand in self.p_ligand:
            model.parameters['Ls_0'].value = ligand

            ligand_array.append(ligand)

            y = odesolve(model, self.tspan, compiler="python")

            if len(observables) > 1:
                # print("{0} + {1}".format(observables[0], observables[1]))
                output_array = y[observables[0]] + y[observables[1]]
                if self.num_samples == 1:
                    np.savetxt("{0}_output".format(observables[0]), y[observables[0]], fmt='%f')
                    np.savetxt("{0}_output".format(observables[1]), y[observables[1]], fmt='%f')
                    np.savetxt("output_array", output_array, fmt='%f')

                output.append(output_array[-1])
            else:
                # print("{0}".format(observables[0]))
                output_array = y[observables[0]]
                output.append(output_array[-1])

        np.savetxt("Ligand_concentrations", ligand_array, fmt='%f')
        np.savetxt("output", output, fmt='%f')


class EarlyPositiveFeedback(PysbTcrSelfWithForeign):
    def __init__(self, steps=8, self_foreign=False, lf=30):
        PysbTcrSelfWithForeign.__init__(self, steps=steps, self_foreign=self_foreign, lf=lf)

    def cycle_4(self, i):

        if i == "Ls":
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

    def add_step_8_sos(self, i):

        if i == "Ls":
            model.parameters['Sos_0'].value = 2000
            Parameter('k_sos_on', 0.0002 / 4)
            Parameter('k_sos_off', 0.005)

        previous_product = self.cycle_4(i)
        product = previous_product + "_Sos"
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
            Parameter('k_sos_on_rgdp', 0.0024)
            Parameter('k_sos_off_rgdp', 3.0)

        add_new_monomer(product_rd)

        Rule('{0}_bind'.format(product_rd),
             eval('{0}()'.format(previous_product)) + Ras_GDP() | eval('{0}()'.format(product_rd)),
             k_sos_on_rgdp, k_sos_off_rgdp)

        product_rt = previous_product + "_Ras_GTP"

        if i == "Ls":
            Parameter('k_sos_on_rgtp', 0.0022)
            Parameter('k_sos_off_rgtp', 0.4)

        add_new_monomer(product_rt)

        Rule('{0}_bind'.format(product_rt),
             eval('{0}()'.format(previous_product)) + Ras_GTP() | eval('{0}()'.format(product_rt)),
             k_sos_on_rgtp, k_sos_off_rgtp)

        return product_rd, product_rt

    def add_step_10(self, i):
        previous_product_rd, previous_product_rt = self.add_step_9(i)
        product_rt_rd = previous_product_rt + "_Ras_GDP"

        if i == "Ls":
            Parameter('k_rgdp_on_sos_rgtp', 0.001)
            Parameter('k_rgdp_off_sos_rgtp', 0.1)
            Parameter('k_cat_3', 0.038 * 1.7)

        add_new_monomer(product_rt_rd)

        Rule('{0}_bind'.format(product_rt_rd),
             eval('{0}()'.format(previous_product_rt)) + Ras_GDP() | eval('{0}()'.format(product_rt_rd)),
             k_rgdp_on_sos_rgtp, k_rgdp_off_sos_rgtp)

        Rule('{0}_cat'.format(product_rt_rd),
             eval('{0}()'.format(product_rt_rd)) >> eval('{0}()'.format(previous_product_rt)) + Ras_GTP(),
             k_cat_3)

        product_rd_rd = previous_product_rd + "_Ras_GDP"

        if i == "Ls":
            Parameter('k_rgdp_on_sos_rgdp', 0.0014)
            Parameter('k_rgdp_off_sos_rgdp', 1.0)
            Parameter('k_cat_4', 0.003)

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
            Parameter('k_rgap_on_rgtp', 0.0348)
            Parameter('k_rgap_off_rgtp', 0.2)
            Parameter('k_cat_5', 0.1)

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
        Model()

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
            elif self.steps == 8:
                product = self.add_step_8_sos(i)
            else:
                product = self.add_step_10(i)

            if "O_{0}".format(product) not in observables:
                observables.append("O_{0}".format(product))

        return observables


class PysbStochasticSimulations(PysbTcrSelfWithForeign):
    def __init__(self, ls, steps=8, self_foreign=False, lf=30):
        PysbTcrSelfWithForeign.__init__(self, steps=steps, self_foreign=self_foreign, lf=lf)
        self.run_time = 1000
        self.tspan = np.linspace(0, self.run_time)
        self.ls = ls
        self.runs = 2

    def main(self):
        observables = self.make_model()
        write_columns(observables)
        write_model_attributes(model.rules, "rules")
        write_model_attributes(model.parameters, "parameters")
        write_model_attributes(model.observables, "observables")
        np.savetxt("time", self.tspan, fmt='%f')

        model.parameters['Ls_0'].value = self.ls
        sim = BngSimulator(model, tspan=self.tspan, cleanup=False)
        simulation_results = sim.run(n_runs=self.runs, method='ssa', cleanup=False,
                                     output_dir=os.getcwd())

        # sim = StochKitSimulator(model, tspan=self.tspan)
        # simulation_results = sim.run(n_runs=self.runs)

        if len(observables) > 1:

            ls_results = []
            lf_results = []
            for i in range(self.runs):
                ls_results.append(simulation_results.observables[i][observables[0]])
                lf_results.append(simulation_results.observables[i][observables[1]])

            output = np.mean(ls_results, axis=0) + np.mean(lf_results, axis=0)
        else:
            results = []
            for i in range(self.runs):
                results.append(simulation_results.observables[i][observables[0]])
            output = np.mean(results, axis=0)
            print(results)
            print(output)

        np.savetxt("output", [output[-1]], fmt='%f')


class PysbTcrCubicFbLoop(PysbTcrSelfWithForeign):
    def __init__(self, steps=8, self_foreign=False, lf=30):
        PysbTcrSelfWithForeign.__init__(self, steps=steps, self_foreign=self_foreign, lf=lf)

    def non_specific_step_8(self, i):

        if i == "Ls":
            if self.p_flag:
                Parameter('k_grb_on', self.parameters['k_8_1'])
                Parameter('k_grb_product', self.parameters['k_8_2'])
            else:
                Parameter('k_grb_on', 0.00005)
                Parameter('k_grb_product', 0.01)

        previous_product = self.non_specific_step_7(i)
        intermediate_product = "LATP_Grb"
        product = "final_product"

        if i == "Ls":
            add_new_monomer(intermediate_product)
            Rule("{0}_convert".format(intermediate_product),
                 eval('{0}()'.format(previous_product)) + Grb() | eval('{0}()'.format(intermediate_product)),
                 k_grb_on, k_p_off_R_pmhc)

            add_new_monomer(product)
            Rule("{0}_cat".format(product),
                 eval('{0}()'.format(intermediate_product)) >> LATP() + eval('{0}()'.format(product)),
                 k_grb_product)
            Rule("{0}_uncat".format(product), eval('{0}()'.format(product)) >> Grb(), k_p_off_R_pmhc)

            add_observable(product)

        return product

    def add_positive_feedback_nonspecific(self, i):

        if i == "Ls":
            if self.p_flag:
                Parameter('k_latp_product', self.parameters['k_9_1'])
                Parameter('k_latp_product_grb', self.parameters['k_9_2'])
                Parameter('k_positive_fb', self.parameters['k_9_3'])
            else:
                Parameter('k_latp_product', 0.0003)
                Parameter('k_latp_product_grb', 0.0004)
                Parameter('k_positive_fb', 5.0)

        previous_product = self.non_specific_step_8(i)
        intermediate_product = "LATP_" + previous_product
        intermediate_product_2 = "LATP_" + previous_product + "_Grb"

        if i == "Ls":
            add_new_monomer(intermediate_product)
            Rule("{0}_bind".format(intermediate_product),
                 LATP() + eval('{0}()'.format(previous_product)) | eval('{0}()'.format(intermediate_product)),
                 k_latp_product, k_p_off_R_pmhc)

            # self.add_observable(intermediate_product)
            add_new_monomer(intermediate_product_2)
            Rule("{0}_bind".format(intermediate_product_2),
                 eval('{0}()'.format(intermediate_product)) + Grb() | eval('{0}()'.format(intermediate_product_2)),
                 k_latp_product_grb, k_p_off_R_pmhc)

            Rule("{0}_cat".format(intermediate_product_2),
                 eval('{0}()'.format(intermediate_product_2)) >> eval('{0}()'.format(intermediate_product)) + eval(
                     '{0}()'.format(previous_product)),
                 k_positive_fb)

        return previous_product


class PysbTcrLckFeedbackLoop(PysbTcrSelfWithForeign):
    def __init__(self):
        PysbTcrSelfWithForeign.__init__(self)

    def cycle_4(self, i):

        if i == "Ls":
            Parameter('k_p_on_zap_species', self.rate_constants.k_p_on_zap_species)
            # Parameter('k_p_off_zap_species', self.rate_constants.k_p_off_zap_species)

        previous_product = self.cycle_3(i)
        product = "RP{0}_Zap_P".format(i)
        add_new_monomer(product)

        Rule('{0}_cat'.format(product), eval('{0}()'.format(previous_product)) >> Lck() + eval('{0}()'.format(product)),
             k_p_on_zap_species)
        add_observable(product)

        no_ligand = self.unbind_ligand(i, product)

        if i == "Ls":
            self.dephosphorylate_zap(no_ligand)

        return product

    def add_negative_feedback(self, i):

        if i == "Ls":
            Parameter('k_negative_fb', self.rate_constants.k_negative_loop)
            Parameter('k_lcki', 0.01)

        previous_product = self.cycle_4(i)
        product = "LckI"
        add_new_monomer(product)

        Rule('{0}_negfb'.format(previous_product), eval('{0}()'.format(previous_product)) + Lck() >> eval(
            '{0}()'.format(product)) + eval('{0}()'.format(previous_product)), k_negative_fb)

        if i == "Ls":
            Rule('{0}_activate'.format(product), LckI() >> Lck(), k_lcki)

        return previous_product

    def cycle_6(self, i):

        if i == "Ls":
            Parameter('k_lat_on_species', self.rate_constants.k_lat_on_species)
            Parameter('k_lat_off_species', self.rate_constants.k_lat_off_species)

        previous_product = self.add_negative_feedback(i)
        product = "RP{0}_Zap_P_LAT".format(i)
        add_new_monomer(product)

        Rule("{0}_bind".format(product), eval('{0}()'.format(previous_product)) + LAT() | eval('{0}()'.format(product)),
             k_lat_on_species, k_lat_off_species)
        add_observable(product)

        no_ligand = self.unbind_ligand(i, product)

        if i == "Ls":
            no_zap_p = self.dephosphorylate_zap(no_ligand)
            self.rp_zap_off(no_zap_p)

        return product


# if __name__ == "__main__":
#
#     parser = argparse.ArgumentParser(description="Submitting ode calculations as function of steps",
#                                      formatter_class=argparse.ArgumentDefaultsHelpFormatter)
#     parser.add_argument('--steps', dest='steps', action='store', type=int, default=8, help="number of KP steps.")
#     # parser.add_argument('--ls_lf', action='store_true', default=False, help="flag for submitting Ls_Lf calculations.")
#     parser.add_argument('--lf', dest='lf', action='store', type=int, help="Number of foreign ligands.")
#     parser.add_argument('--ssa', dest='ssa', action='store', type=int, help="Run SSA simulations. Specify amount of Ls.")
#     parser.add_argument('--early_pos_fb', dest='early_pos_fb', action='store_true', default=False,
#                         help='Flag for building and submitting early positive feedback loop.')
#
#     args = parser.parse_args()
#
#     if args.lf:
#         if args.early_pos_fb:
#             tcr = EarlyPositiveFeedback(steps=args.steps, self_foreign=True, lf=args.lf)
#
#         elif args.ssa:
#             tcr = PysbStochasticSimulations(args.ssa, steps=args.steps, self_foreign=True, lf=args.lf)
#
#         else:
#             tcr = PysbTcrSelfWithForeign(steps=args.steps, self_foreign=True, lf=args.lf)
#
#         # if args.cubic:
#         #     tcr = PysbTcrCubicFbLoop(steps=args.steps, self_foreign=True, lf=args.lf)
#     else:
#         if args.early_pos_fb:
#             tcr = EarlyPositiveFeedback(steps=args.steps)
#
#         elif args.ssa:
#             tcr = PysbStochasticSimulations(args.ssa, steps=args.steps)
#
#         else:
#             tcr = PysbTcrSelfWithForeign(steps=args.steps)
#         # if args.cubic:
#         #     tcr = PysbTcrCubicFbLoop(steps=args.steps)
#
#     tcr.main()

# Uncomment to make reaction network

tcr = PysbTcrSelfWithForeign(steps=7, self_foreign=True)
tcr.make_model()

# tcr = PysbTcrCubicFbLoop(steps=8, self_foreign=True)
# tcr.make_model()
# tcr = PysbTcrCubicFbLoop(steps=8, self_foreign=True)
# tcr.make_model()

### python /anaconda2/lib/python2.7/site-packages/pysb/tools/render_reactions.py ../SSC_python_modules/pysb_t_cell_network.py > model.dot
### dot -Tpdf model.dot -o model.pdf
