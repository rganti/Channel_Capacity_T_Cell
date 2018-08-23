import argparse

import numpy as np
from pysb import *
from pysb.integrate import odesolve

from realistic_network import BindingParameters


def e(i, s=""):
    eval('{0}{1}()'.format(i, s))


class PysbTcrSelfWithForeign(object):
    def __init__(self, steps=8, self_foreign=False, lf=30):
        self.rate_constants = BindingParameters()

        self.lf = lf
        self.ls = 1000
        self.steps = steps

        self.run_time = 10000
        self.tspan = np.linspace(0, self.run_time)
        self.simulation_name = "ODE_" + str(self.steps)

        self.mu = 6
        self.sigma = 1.0
        self.num_samples = 1000

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

        Parameter('R_0', 10000)
        Parameter('Ls_0', self.ls)
        Parameter('Lck_0', 2000)
        Parameter('Zap_0', 3000)
        Parameter('LAT_0', 1000)

        Initial(R(), R_0)
        Initial(Ls(), Ls_0)
        Initial(Lck(), Lck_0)
        Initial(Zap(), Zap_0)
        Initial(LAT(), LAT_0)

        if "Lf" in self.output:
            Monomer('Lf')
            Parameter('Lf_0', self.lf)
            Initial(Lf(), Lf_0)

    def define_rates(self):
        Parameter('k_L_on', self.rate_constants.k_L_on)
        if "Lf" in self.output:
            Parameter('k_Lf_off', self.rate_constants.k_foreign_off)
        Parameter('k_Ls_off', self.rate_constants.k_self_off)

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

    def rp_zap_off(self, product):
        new = product.replace("RP_Zap_", "")
        self.add_new_monomer(new)
        Rule('{0}_rpzap_off'.format(product), eval('{0}()'.format(product)) | RP_Zap() + eval('{0}()'.format(new)),
             k_zap_off_R, k_zap_on_R)

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
            Parameter('k_lat_off_species', self.rate_constants.k_lat_off_species)

        previous_product = self.add_negative_feedback(i)
        product = "RP{0}_Lck_Zap_P_LAT".format(i)
        self.add_new_monomer(product)

        Rule("{0}_bind".format(product), eval('{0}()'.format(previous_product)) + LAT() | eval('{0}()'.format(product)),
             k_lat_on_species, k_lat_off_species)
        self.add_observable(product)

        no_ligand = self.unbind_ligand(i, product)

        if i == "Ls":
            no_lck = self.lck_off(no_ligand, k_off="k_lck_off_zap_R")
            no_zap_p = self.dephosphorylate_zap(no_lck)
            self.rp_zap_off(no_zap_p)

        return product

    def non_specific_step_7(self, i):
        if i == "Ls":
            Parameter('k_p_lat_on_species', self.rate_constants.k_p_lat_on_species)

        previous_product = self.cycle_6(i)
        product = "LATP"
        if i == "Ls":
            self.add_new_monomer(product)

        Rule("{0}_cat".format(previous_product),
             eval('{0}()'.format(previous_product)) >> eval('{0}()'.format("RP{0}_Lck_Zap_P".format(i))) + eval(
                 '{0}()'.format(product)), k_p_lat_on_species)
        if i == "Ls":
            self.add_observable(product)

            Rule("{0}_uncat".format(product), eval('{0}()'.format(product)) >> LAT(), k_p_off_R_pmhc)

        return product

    def add_step_7(self, i):

        if i == "Ls":
            Parameter('k_p_lat_on_species', self.rate_constants.k_p_lat_on_species)

        previous_product = self.cycle_6(i)
        product = "{0}_LATP".format(i)
        self.add_new_monomer(product)

        Rule("{0}_cat".format(product),
             eval('{0}()'.format(previous_product)) >> eval('{0}()'.format("RP{0}_Lck_Zap_P".format(i))) + eval(
                 '{0}()'.format(product)),
             k_p_lat_on_species)
        self.add_observable(product)

        Rule("{0}_uncat".format(product), eval('{0}()'.format(product)) >> LAT(), k_p_off_R_pmhc)

        return product

    def non_specific_step_8(self, i):

        if i == "Ls":
            Parameter('k_product_on', self.rate_constants.k_product_on)

        previous_product = self.non_specific_step_7(i)
        product = "final_product"
        if i == "Ls":
            self.add_new_monomer(product)
            Rule("{0}_convert".format(product), eval('{0}()'.format(previous_product)) | eval('{0}()'.format(product)),
                 k_product_on, k_p_off_R_pmhc)
            self.add_observable(product)

        return product

    def add_step_8(self, i):

        if i == "Ls":
            Parameter('k_product_on', self.rate_constants.k_product_on)

        previous_product = self.add_step_7(i)
        product = "{0}_product".format(i)
        self.add_new_monomer(product)

        Rule("{0}_convert".format(product), eval('{0}()'.format(previous_product)) | eval('{0}()'.format(product)),
             k_product_on, k_p_off_R_pmhc)
        self.add_observable(product)

        return product

    def add_positive_feedback_nonspecific(self, i):

        if i == "Ls":
            Parameter('k_positive_fb', self.rate_constants.k_positive_loop)

        product = self.non_specific_step_8(i)
        if i == "Ls":
            Rule("{0}_posfb".format(product),
                 eval('{0}()'.format(product)) + eval('{0}()'.format("LATP")) >> eval(
                     '{0}()'.format(product)) + eval('{0}()'.format(product)),
                 k_positive_fb)

        return product

    def add_positive_feedback(self, i):

        if i == "Ls":
            Parameter('k_positive_fb', self.rate_constants.k_positive_loop)

        product = self.add_step_8(i)

        Rule("{0}_posfb".format(product),
             eval('{0}()'.format(product)) + eval('{0}()'.format("{0}_LATP".format(i))) >> eval(
                 '{0}()'.format(product)) + eval('{0}()'.format(product)),
             k_positive_fb)

        return product

    def write_columns(self, observables):
        f = open("column_names", "w")
        for item in observables:
            f.write("{0} ".format(item))
        f.write("\n")
        f.close()

    def write_model_attributes(self, attributes, file_name):
        f = open(file_name, "w")
        for item in attributes:
            f.write("{0} \n".format(item))
        f.write("\n")
        f.close()

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
                # product = self.add_step_7(i)
                product = self.non_specific_step_7(i)
            elif self.steps == 8:
                # product = self.add_step_8(i)
                product = self.non_specific_step_8(i)
            else:
                # product = self.add_positive_feedback(i)
                product = self.add_positive_feedback_nonspecific(i)

            if "O_{0}".format(product) not in observables:
                observables.append("O_{0}".format(product))

        return observables

    def main(self):
        output = []
        ligand_array = []
        observables = self.make_model()
        # print(observables)

        self.write_columns(observables)
        self.write_model_attributes(model.rules, "rules")
        self.write_model_attributes(model.parameters, "parameters")
        self.write_model_attributes(model.observables, "observables")

        np.savetxt("time", self.tspan, fmt='%f')

        for ligand in self.p_ligand:
            model.parameters['Ls_0'].value = ligand

            ligand_array.append(ligand)

            y = odesolve(model, self.tspan)

            if len(observables) > 1:
                print("{0} + {1}".format(observables[0], observables[1]))
                output_array = y[observables[0]] + y[observables[1]]
                if self.num_samples == 1:
                    np.savetxt("{0}_output".format(observables[0]), y[observables[0]], fmt='%f')
                    np.savetxt("{0}_output".format(observables[1]), y[observables[1]], fmt='%f')
                    np.savetxt("output_array", output_array, fmt='%f')

                output.append(output_array[-1])
            else:
                print("{0}".format(observables[0]))
                output_array = y[observables[0]]
                output.append(output_array[-1])

        np.savetxt("Ligand_concentrations", ligand_array, fmt='%f')
        np.savetxt("output", output, fmt='%f')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Submitting ode calculations as function of steps",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--steps', dest='steps', action='store', type=int, default=8,
                        help="number of KP steps.")
    parser.add_argument('--ls_lf', action='store_true', default=False,
                        help="flag for submitting Ls_Lf calculations.")
    parser.add_argument('--lf', dest='lf', action='store', type=int, default=30,
                        help="number of foreign ligands.")

    args = parser.parse_args()

    if args.ls_lf:
        tcr = PysbTcrSelfWithForeign(steps=args.steps, self_foreign=True, lf=args.lf)
    else:
        tcr = PysbTcrSelfWithForeign(steps=args.steps)

    tcr.main()

# tcr = PysbTcrSelfWithForeign(steps=8, self_foreign=True)
# tcr.make_model()
