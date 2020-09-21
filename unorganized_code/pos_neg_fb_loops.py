import os

import numpy as np
from pysb import Parameter, Rule
from pysb.simulator import BngSimulator

from pysb_t_cell_network import PysbTcrSelfWithForeign, add_observable, write_columns, \
    write_model_attributes


class EarlyPositiveFeedback(PysbTcrSelfWithForeign):
    def __init__(self, steps=8, self_foreign=False, lf=30):
        PysbTcrSelfWithForeign.__init__(self, steps=steps, self_foreign=self_foreign, lf=lf)

    def cycle_4(self, i):

        if i == "Ls":
            Parameter('k_p_on_zap_species', self.rate_constants.k_p_on_zap_species)
            Parameter('k_p_off_zap_species', self.rate_constants.k_p_off_zap_species)

        previous_product = self.cycle_3(i)
        product = "RP{0}_Lck_Zap_P".format(i)
        self.add_new_monomer(product)

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
            self.model.parameters['Sos_0'].value = 2000
            Parameter('k_sos_on', 0.0002 / 4)
            Parameter('k_sos_off', 0.005)

        previous_product = self.cycle_4(i)
        product = previous_product + "_Sos"
        self.add_new_monomer(product)

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

        self.add_new_monomer(product_rd)

        Rule('{0}_bind'.format(product_rd),
             eval('{0}()'.format(previous_product)) + Ras_GDP() | eval('{0}()'.format(product_rd)),
             k_sos_on_rgdp, k_sos_off_rgdp)

        product_rt = previous_product + "_Ras_GTP"

        if i == "Ls":
            Parameter('k_sos_on_rgtp', 0.0022)
            Parameter('k_sos_off_rgtp', 0.4)

        self.add_new_monomer(product_rt)

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

        self.add_new_monomer(product_rt_rd)

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

        self.add_new_monomer(product_rd_rd)

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

            self.add_new_monomer(product)

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
        write_model_attributes(self.model.rules, "rules")
        write_model_attributes(self.model.parameters, "parameters")
        write_model_attributes(self.model.observables, "observables")
        np.savetxt("time", self.tspan, fmt='%f')

        self.model.parameters['Ls_0'].value = self.ls
        sim = BngSimulator(self.model, tspan=self.tspan, cleanup=False)
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

        previous_product = self.step_8(i)
        intermediate_product = "LATP_Grb"
        product = "final_product"

        if i == "Ls":
            self.add_new_monomer(intermediate_product)
            Rule("{0}_convert".format(intermediate_product),
                 eval('{0}()'.format(previous_product)) + Grb() | eval('{0}()'.format(intermediate_product)),
                 k_grb_on, k_p_off_R_pmhc)

            self.add_new_monomer(product)
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
            self.add_new_monomer(intermediate_product)
            Rule("{0}_bind".format(intermediate_product),
                 LATP() + eval('{0}()'.format(previous_product)) | eval('{0}()'.format(intermediate_product)),
                 k_latp_product, k_p_off_R_pmhc)

            # self.add_observable(intermediate_product)
            self.add_new_monomer(intermediate_product_2)
            Rule("{0}_bind".format(intermediate_product_2),
                 eval('{0}()'.format(intermediate_product)) + Grb() | eval('{0}()'.format(intermediate_product_2)),
                 k_latp_product_grb, k_p_off_R_pmhc)

            Rule("{0}_cat".format(intermediate_product_2),
                 eval('{0}()'.format(intermediate_product_2)) >> eval('{0}()'.format(intermediate_product)) + eval(
                     '{0}()'.format(previous_product)),
                 k_positive_fb)

        return previous_product
