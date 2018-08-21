# import datetime

# import matplotlib.pyplot as plt
import numpy as np
from pysb import *
from pysb.integrate import odesolve

from realistic_network import BindingParameters


def e(i, s=""):
    eval('{0}{1}()'.format(i, s))


class PysbTcrSelfWithForeign(object):
    def __init__(self):
        self.rate_constants = BindingParameters()

        self.lf = 30
        self.ls = 800

        self.run_time = 100
        self.tspan = np.linspace(0, self.run_time)
        self.simulation_time = 1
        self.simulation_name = "ODE_Neg_FB"

        self.output = ["Ls", "Lf"]
        # self.output = ["Ls"]

    def define_monomers(self):
        Monomer('R')
        Monomer('Ls')
        Monomer('Lck')
        Monomer('Zap')

        Parameter('R_0', 10000)
        Parameter('Ls_0', self.ls)
        Parameter('Lck_0', 2000)
        Parameter('Zap_0', 3000)

        Initial(R(), R_0)
        Initial(Ls(), Ls_0)
        Initial(Lck(), Lck_0)
        Initial(Zap(), Zap_0)

        if "Lf" in self.output:
            Monomer('Lf')
            Parameter('Lf_0', self.lf)
            Initial(Lf(), Lf_0)

    def define_rates(self):
        Parameter('k_L_on', self.rate_constants.k_L_on)
        if "Lf" in self.output:
            Parameter('k_Lf_off', self.rate_constants.k_foreign_off)
        Parameter('k_Ls_off', self.rate_constants.k_self_off)

    def ligand_specificity(self, i):
        if i == "Ls":
            l_unbound = Ls(b=None)
            l_bound = Ls(b=1)
        else:
            l_unbound = Lf(b=None)
            l_bound = Lf(b=1)

        return l_unbound, l_bound

    def add_new_monomer(self, product):
        try:
            model.monomers[product]
        except:
            Monomer(product)
            Initial(eval(product + "()"), Parameter(product + "_0", 0))

    def add_observable(self, species):
        Observable("O_{0}".format(species), eval('{0}()'.format(species)))

    def add_step_0(self, i):
        product = "R{0}".format(i)
        self.add_new_monomer(product)

        Rule('{0}_bind'.format(product), R() + eval('{0}()'.format(i)) | eval('{0}()'.format(product)),
             k_L_on, eval('k_{0}_off'.format(i)))

        self.add_observable(i)
        self.add_observable(product)
        return product

        # if i == "Ls":
        #     # Initial(RLs(), Parameter('RLs_0', 0))
        #     Rule('R_{0}_bind'.format(i), R() + Ls() | RLs(), k_L_on, k_Ls_off)
        #     Observable("O_{0}".format(i), Ls())
        #     Observable("O_R{0}".format(i), RLs())
        # else:
        #     Initial(RLf(), Parameter('RLf_0', 0))
        #     Rule('R_{0}_bind'.format(i), R() + Lf() | RLf(), k_L_on, k_Lf_off)
        #     Observable("O_{0}".format(i), Lf())
        #     Observable("O_R{0}".format(i), RLf())
        #
        # Observable("O_{0}".format(i), l_unbound)
        # Observable("O_R{0}".format(i), R(b=1, lck=None) % l_bound)

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
            Parameter('k_lck_on_RL', self.rate_constants.k_lck_on_R_pmhc)
            Parameter('k_lck_off_RL', self.rate_constants.k_lck_off_R_pmhc)

            Parameter('k_lck_on_R', self.rate_constants.k_lck_on_R)
            Parameter('k_lck_off_R', self.rate_constants.k_lck_off_R)

        previous_product = self.add_step_0(i)

        product = "R{0}_Lck".format(i)
        self.add_new_monomer(product)
        Rule('{0}_bind'.format(product), eval('{0}()'.format(previous_product)) + Lck() | eval('{0}()'.format(product)),
             k_lck_on_RL, k_lck_off_RL)
        self.add_observable(product)

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
        self.add_new_monomer(product)
        Rule('{0}_cat'.format(product), eval('{0}()'.format(previous_product)) | eval('{0}()'.format(product)),
             k_p_on_R_pmhc, k_p_off_R_pmhc)
        self.add_observable(product)

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
        self.add_new_monomer(product)

        Rule('{0}_bind'.format(product), eval('{0}()'.format(previous_product)) + Zap() | eval('{0}()'.format(product)),
             k_zap_on_R_pmhc, k_zap_off_R_pmhc)
        self.add_observable(product)

        no_ligand = self.unbind_ligand(i, product)

        if i == "Ls":
            no_lck = self.lck_off(no_ligand, k_off="k_lck_off_zap_R")
            no_zap = self.zap_off(no_lck)

        return product

    def dephosphorylate_zap(self, product):
        new = product.replace("_Zap_P", "_Zap")
        self.add_new_monomer(new)
        Rule('{0}_uncat'.format(product), eval('{0}()'.format(product)) | eval('{0}()'.format(new)),
             k_p_off_R, k_p_on_R)

        return new

    def cycle_4(self, i):

        if i == "Ls":
            Parameter('k_p_on_zap_species', self.rate_constants.k_p_on_zap_species)
            Parameter('k_p_off_zap_species', self.rate_constants.k_p_off_zap_species)

        previous_product = self.cycle_3(i)
        product = "RP{0}_Lck_Zap_P".format(i)
        self.add_new_monomer(product)

        Rule('{0}_cat'.format(product), eval('{0}()'.format(previous_product)) | eval('{0}()'.format(product)),
             k_p_on_zap_species, k_p_off_zap_species)
        self.add_observable(product)

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
        self.add_new_monomer("R{0}_LckI".format(i))

        Rule('{0}_negfb_1'.format(previous_product),
             eval('{0}()'.format(previous_product)) + eval('R{0}_Lck()'.format(i)) >> eval(
                 'R{0}_LckI()'.format(i)) + eval('{0}()'.format(previous_product)),
             k_negative_fb)

        self.add_new_monomer("LckI")

        Rule('{0}_negfb_2'.format('R{0}_LckI'.format(i)),
             eval('R{0}_LckI()'.format(i)) >> eval('R{0}()'.format(i)) + eval('LckI()'),
             k_diss_rllcki)

        if i == "Ls":
            Rule('{0}_negfb_3'.format("LckI"), LckI() >> Lck(), k_lcki)

        return previous_product

    def cycle_6(self, i):

        if i == "Ls":
            Parameter('k_lat_on_species', self.rate_constants.k_lat_on_species)

        # # Sixth LAT on
        self.k_lat_on_species = self.k_lck_on_R_pmhc / 100.0
        self.k_lat_off_species = self.k_lck_off_R_pmhc

        self.k_lat_on_rp_zap = self.k_zap_on_R
        self.k_lat_off_rp_zap = 20.0

    def main(self):

        Model()

        self.define_monomers()
        self.define_rates()

        x = []
        for i in self.output:
            # product = self.add_step_0(i)
            # product = self.cycle_1(i)
            # product = self.cycle_2(i)
            # product = self.cycle_3(i)
            # product = self.cycle_4(i)
            product = self.add_negative_feedback(i)
            x.append("O_{0}".format(product))

        print(x)

        y = odesolve(model, self.tspan)

        np.savetxt("output", np.array([y[x[0]], y[x[1]]]), fmt='%f')

        #
        # plt.plot(self.tspan, y['O_RLf'], label="RLf")
        # plt.plot(self.tspan, y['O_RLs'], label="RLs")
        # plt.legend()
        # plt.savefig("ode.pdf", format='pdf')


tcr = PysbTcrSelfWithForeign()
tcr.main()
