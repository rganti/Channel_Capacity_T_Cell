import argparse
import datetime
import os
import subprocess

import numpy as np
import pandas as pd
from pysb import *
from pysb.integrate import odesolve

from pysb_t_cell_network import write_columns, write_model_attributes
from realistic_network import make_and_cd


def add_new_monomer(product):
    try:
        model.monomers[product]
    except:
        Monomer(product)
        Initial(eval(product + "()"), Parameter(product + "_0", 0))


def add_observable(species):
    Observable("O_{0}".format(species), eval('{0}()'.format(species)))


class SoSFeedback(object):

    def __init__(self):
        self.run_time = 300
        self.tspan = np.linspace(0, self.run_time)

        self.sos = [round(i) for i in np.linspace(25, 500, num=40)]

        self.model = Model()

    def define_monomers(self):
        Monomer('Sos')
        Monomer('Ras_GDP')
        Monomer('Ras_GTP')
        Monomer('Ras_GAP')

        Parameter('Sos_0', self.sos[0])
        Parameter('Ras_GDP_0', 300)
        Parameter('Ras_GAP_0', 10)
        Parameter('Ras_GTP_0', 0)

        Initial(Sos(), Sos_0)
        Initial(Ras_GDP(), Ras_GDP_0)
        Initial(Ras_GAP(), Ras_GAP_0)
        Initial(Ras_GTP(), Ras_GTP_0)

    def add_step_1(self):
        Parameter('k_sos_on_rgdp', 0.0024)
        Parameter('k_sos_off_rgdp', 3.0)

        product = "Sos_Ras_GDP"
        add_new_monomer(product)

        Rule('{0}_bind'.format(product), Sos() + Ras_GDP() | eval('{0}()'.format(product)),
             k_sos_on_rgdp, k_sos_off_rgdp)

        add_observable(product)
        return product

    def add_step_2(self):
        Parameter('k_sos_on_rgtp', 0.0022)
        Parameter('k_sos_off_rgtp', 0.4)

        product = "Sos_Ras_GTP"
        add_new_monomer(product)

        Rule('{0}_bind'.format(product), Sos() + Ras_GTP() | eval('{0}()'.format(product)),
             k_sos_on_rgtp, k_sos_off_rgtp)

        add_observable(product)

        return product

    def add_step_3(self):
        Parameter('k_rgdp_on_sos_rgtp', 0.001)
        Parameter('k_rgdp_off_sos_rgtp', 0.1)
        Parameter('k_cat_3', 0.038 * 1.7)

        previous_product = self.add_step_2()

        product = "Sos_Ras_GTP_Ras_GDP"
        add_new_monomer(product)

        Rule('{0}_bind'.format(product),
             eval('{0}()'.format(previous_product)) + Ras_GDP() | eval('{0}()'.format(product)),
             k_rgdp_on_sos_rgtp, k_rgdp_off_sos_rgtp)

        Rule('{0}_cat'.format(product),
             eval('{0}()'.format(product)) >> eval('{0}()'.format(previous_product)) + Ras_GTP(),
             k_cat_3)

    def add_step_4(self):
        Parameter('k_rgdp_on_sos_rgdp', 0.0014)
        Parameter('k_rgdp_off_sos_rgdp', 1.0)
        Parameter('k_cat_4', 0.003)

        previous_product = self.add_step_1()

        product = "Sos_Ras_GDP_Ras_GDP"
        add_new_monomer(product)

        Rule('{0}_bind'.format(product),
             eval('{0}()'.format(previous_product)) + Ras_GDP() | eval('{0}()'.format(product)),
             k_rgdp_on_sos_rgdp, k_rgdp_off_sos_rgdp)

        Rule('{0}_cat'.format(product),
             eval('{0}()'.format(product)) >> eval('{0}()'.format(previous_product)) + Ras_GTP(),
             k_cat_4)

    def add_step_5(self):
        Parameter('k_rgap_on_rgtp', 0.0348)
        Parameter('k_rgap_off_rgtp', 0.2)
        Parameter('k_cat_5', 0.1)

        product = "Ras_GAP_Ras_GTP"
        add_new_monomer(product)

        Rule('{0}_bind'.format(product),
             Ras_GAP() + Ras_GTP() | eval('{0}()'.format(product)),
             k_rgap_on_rgtp, k_rgap_off_rgtp)

        add_observable(product)

        Rule('{0}_cat'.format(product),
             eval('{0}()'.format(product)) >> Ras_GAP() + Ras_GDP(),
             k_cat_5)

    def make_model(self):
        # Model()

        self.define_monomers()

        observables = []

        self.add_step_3()
        self.add_step_4()
        self.add_step_5()

        add_observable("Ras_GTP")

        product = "Ras_GTP"
        observables.append("O_{0}".format(product))

        return observables

    def main(self):
        sos_array = []
        output = []
        observables = self.make_model()

        write_columns(observables)
        write_model_attributes(model.rules, "rules")
        write_model_attributes(model.parameters, "parameters")
        write_model_attributes(model.observables, "observables")

        np.savetxt("time", self.tspan, fmt='%f')

        for sos in self.sos:
            model.parameters['Sos_0'].value = sos
            y = odesolve(model, self.tspan, compiler="python")

            sos_array.append(sos)
            # print(y[observables[0]][-1])
            output.append(y[observables[0]][-1])

        df = pd.DataFrame({'Sos': sos_array, 'RasGTP': output})
        df.to_csv("./sos_rasgtp", sep='\t')

        # np.savetxt("Sos", sos_array, fmt='%f')
        # np.savetxt("RasGTP", output, fmt='%f')


class SoSFeedbackLigandSpecific(SoSFeedback):
    def __init__(self):
        SoSFeedback.__init__(self)

        self.sos_total = [round(i) for i in np.linspace(25, 500, num=40)]


class LaunchQsub(object):
    def __init__(self):
        self.simulation_name = "Sos_FB"
        self.simulation_time = 2

    def generate_qsub(self):
        q = open("qsub.sh", "w")
        q.write("#PBS -m ae\n")
        q.write("#PBS -q short\n")
        q.write("#PBS -V\n")
        q.write("#PBS -l walltime={1},nodes=1:ppn=2 -N {0}\n\n".format(self.simulation_name,
                                                                       datetime.timedelta(
                                                                           minutes=self.simulation_time)))
        q.write("cd $PBS_O_WORKDIR\n")
        q.write("echo $PBS_JOBID > job_id\n\n")

        q.write(
            "python ~/SSC_python_modules/pysb_sos_fb.py --run\n")
        q.close()

    def launch(self):
        (stdout, stderr) = subprocess.Popen(["qsub {0}".format("qsub.sh")], shell=True, stdout=subprocess.PIPE,
                                            cwd=os.getcwd()).communicate()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Submitting ode calculations as function of steps",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--run', action='store_true', default=False,
                        help='Flag for submitting simulations.')

    args = parser.parse_args()

    if not args.run:
        make_and_cd("Ras_SoS_Fb_200")
        qsub = LaunchQsub()
        qsub.generate_qsub()
        qsub.launch()
    else:
        sos = SoSFeedback()
        sos.main()
