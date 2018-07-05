#!/usr/bin/python
from two_species import SelfWithForeign, compile_script


class SimpleSecondOrder(SelfWithForeign):
    def __init__(self):
        SelfWithForeign.__init__(self)

        self.n_initial = {"R": 1000, "Lf": 500}
        self.record = ["Lf", "C0"]

        self.simulation_name = "second_order"

        self.forward_rates = {"RLf": self.rate_constants.kappa}
        self.reverse_rates = {"C0": self.rate_constants.k_off_foreign}

        self.forward_rxns = [[["R", "Lf"], ["C0"]]]
        self.reverse_rxns = [[["C0"], ["R", "Lf"]]]


if __name__ == "__main__":

    second_order = SimpleSecondOrder()

    second_order.generate_ssc_script()
    compile_script(second_order.simulation_name + ".rxn")
    second_order.generate_qsub()










