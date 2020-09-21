from pysb import *

from pysb_t_cell_network import PysbTcrSelfWithForeign


class PysbTcrSelfWithForeignJupyter(PysbTcrSelfWithForeign):

    def __init__(self):
        PysbTcrSelfWithForeign.__init__(self)

    def main(self):
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
                product = self.add_step_7(i)
            elif self.steps == 8:
                product = self.add_step_8(i)
            else:
                product = self.add_positive_feedback(i)

            observables.append("O_{0}".format(product))


tcr = PysbTcrSelfWithForeignJupyter()

tcr.main()
