import argparse

import numpy as np
from pysb import *
from pysb.integrate import odesolve

from pysb_t_cell_network import write_columns, write_model_attributes


class NewBindingParameters(object):
    def __init__(self):
        # Zeroth Cycle ligand binding
        self.k_L_on = 0.0022
        self.k_foreign_off = 0.2
        self.k_self_off = 2.0
        self.kp = 1.0
        self.k_off = 0.5


class ToyModel(object):
    def __init__(self, self_foreign=False):
        self.rate_constants = NewBindingParameters()

        self.steps = 3
        self.model = Model()
        self.self_foreign = self_foreign

        if self.self_foreign:
            self.output = ["Ls", "Lf"]
        else:
            self.output = ["Ls"]

        self.mu = 5.0
        self.sigma = 1.0
        self.num_samples = 1000

        self.run_time = 10000
        self.tspan = np.linspace(0, self.run_time)

        self.p_ligand = [int(i) for i in np.round(np.random.lognormal(self.mu, self.sigma, self.num_samples))]

    def define_monomers(self):
        Monomer('R')
        Monomer('Ls')
        Monomer('Lf')

        Parameter('R_0', 10000)
        self.add_observable('R')

        Parameter('Ls_0', 30)
        Parameter('Lf_0', 30)

        Initial(R(), R_0)
        Initial(Ls(), Ls_0)
        Initial(Lf(), Lf_0)

    def add_new_monomer(self, product):
        try:
            self.model.monomers[product]
        except:
            Monomer(product)
            Initial(eval(product + "()"), Parameter(product + "_0", 0))

    def add_observable(self, species):
        Observable("O_{0}".format(species), eval('{0}()'.format(species)))

    def define_rates(self):
        Parameter('k_L_on', self.rate_constants.k_L_on)
        Parameter('k_Lf_off', self.rate_constants.k_foreign_off)
        Parameter('k_Ls_off', self.rate_constants.k_self_off)
        Parameter('kp', self.rate_constants.kp)
        Parameter('k_off', self.rate_constants.k_off)

    def add_step_0(self, i):
        product = "R{0}".format(i)
        self.add_new_monomer(product)

        Rule('{0}_bind'.format(product), R() + eval('{0}()'.format(i)) | eval('{0}()'.format(product)),
             k_L_on, eval('k_{0}_off'.format(i)))

        self.add_observable(i)
        self.add_observable(product)
        return product

    def add_cycle(self, i, previous_product, product):
        self.add_new_monomer(product)

        Rule('{0}_bind'.format(product), eval('{0}()'.format(previous_product)) >> eval('{0}()'.format(product)),
             kp)
        self.add_observable(product)

        Rule('{0}_return'.format(product), eval('{0}()'.format(product)) >> R() + eval('{0}()'.format(i)),
             eval('k_{0}_off'.format(i)))

    def cycle_1(self, i):
        previous_product = self.add_step_0(i)
        product = "RP1{0}".format(i)
        self.add_cycle(i, previous_product, product)

        return product

    def cycle_2(self, i):
        previous_product = self.cycle_1(i)
        product = "RP2{0}".format(i)
        self.add_cycle(i, previous_product, product)

        return product

    def cycle_3(self, i):
        previous_product = self.cycle_2(i)
        product = "RP3"

        self.add_new_monomer(product)

        Rule('{0}_{1}_bind'.format(product, i),
             eval('{0}()'.format(previous_product)) >> eval('{0}()'.format(product)) + eval('{0}()'.format(i)), kp)

        if i == "Ls":
            self.add_observable(product)
            Rule('{0}_return'.format(product), eval('{0}()'.format(product)) >> R(),
                 eval('k_off'))

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
            else:
                product = self.cycle_3(i)

            if "O_{0}".format(product) not in observables:
                observables.append("O_{0}".format(product))

        return observables

    def main(self):
        output_array = []
        ligand_array = []
        ls_ss_array = []
        r_ss_array = []
        lf_ss_array = []

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

            output = y[observables[0]]
            ls_ss = y['O_Ls']
            r_ss = y['O_R']

            if self.self_foreign:
                lf_ss = y['O_Lf']
                lf_ss_array.append(lf_ss[-1])

            output_array.append(output[-1])
            ls_ss_array.append(ls_ss[-1])
            r_ss_array.append(r_ss[-1])

        np.savetxt("Ligand_concentrations", ligand_array, fmt='%f')
        np.savetxt("output", output_array, fmt='%f')
        np.savetxt("ls_ss", ls_ss_array, fmt='%f')
        np.savetxt("lf_ss", lf_ss_array, fmt='%f')
        np.savetxt("r_ss", r_ss_array, fmt='%f')


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Submitting ode calculations as function of steps",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--lf', dest='lf', action='store_true', default=False,
                        help="Flag to submit self w/ foreign sims.")
    args = parser.parse_args()

    if args.lf:
        tcr = ToyModel(self_foreign=True)

    else:
        tcr = ToyModel()

    tcr.main()
