import argparse
import datetime
import os
import pickle
import subprocess

import numpy as np
import pandas as pd
from sklearn.model_selection import ParameterGrid

from realistic_network import make_and_cd


class ParameterTesting(object):
    def __init__(self, steps=8):

        self.k_7 = [0.01]

        self.k_8_1 = np.linspace(5.0e-5, 1.0e-4, 5)
        self.k_8_2 = [np.linspace(0.001, 0.01, 3)[-1]]

        self.k_9_1 = np.linspace(0.00005, 0.0003, 3)
        self.k_9_2 = np.linspace(0.00005, 0.0003, 3)
        self.k_9_3 = [100 * i for i in self.k_8_2]

        self.steps = steps

        if self.steps < 9:
            param_grid = {'k_7': [0.01], 'k_8_1': np.linspace(5.0e-5, 1.0e-4, 5),
                          'k_8_2': [np.linspace(0.001, 0.01, 3)[-1]]}

        else:
            param_grid = {'k_7': [0.01],
                          'k_8_1': np.linspace(5.0e-5, 1.0e-4, 5),
                          'k_8_2': [np.linspace(0.001, 0.01, 3)[-1]],
                          'k_9_1': np.linspace(0.00005, 0.0003, 3),
                          'k_9_2': np.linspace(0.00005, 0.0003, 3),
                          'k_9_3': [100 * i for i in self.k_8_2]}

        self.grid = ParameterGrid(param_grid)

        # self.tcr_self_foreign = PysbTcrCubicFbLoop(steps=self.steps, self_foreign=True)
        # self.tcr_self = PysbTcrCubicFbLoop(steps=self.steps)
        self.simulation_name = "ODE_steps_" + str(self.steps)
        self.simulation_time = 10
        self.self_foreign = True

    def generate_qsub(self):
        q = open("qsub.sh", "w")
        q.write("#PBS -m ae\n")
        q.write("#PBS -q short\n")
        q.write("#PBS -V\n")
        q.write("#PBS -l walltime={1},nodes=1:ppn=4 -N {0}\n\n".format(self.simulation_name,
                                                                       datetime.timedelta(
                                                                           minutes=self.simulation_time)))
        q.write("cd $PBS_O_WORKDIR\n\n")
        q.write("echo $PBS_JOBID > job_id\n")

        if self.self_foreign:
            q.write(
                "python ~/SSC_python_modules/pysb_t_cell_network.py --steps {0} --ls_lf --lf {1}\n".format(self.steps,
                                                                                                           self.lf))
        else:
            q.write("python ~/SSC_python_modules/pysb_t_cell_network.py --steps {0}".format(self.steps))
        q.close()

    def launch(self):
        (stdout, stderr) = subprocess.Popen(["qsub {0}".format("qsub.sh")], shell=True, stdout=subprocess.PIPE,
                                            cwd=os.getcwd()).communicate()

    def run_tests(self):
        d = {'file_path': []}
        count = 0
        home_directory = os.getcwd()
        for set in list(self.grid):
            print(str(set))
            file_path = "{0}_step_{1}".format(self.steps, count)
            d['file_path'].append(file_path)
            make_and_cd(file_path)
            pickle_out = open("parameters.pickle", "wb")
            pickle.dump(set, pickle_out)
            pickle_out.close()
            os.chdir(home_directory)

            count += 1

        df = pd.DataFrame(data=d)
        df.to_csv("./file_paths", sep='\t')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Submitting ode calculations as function of steps",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--steps', dest='steps', action='store', type=int, default=8,
                        help="number of KP steps.")

    args = parser.parse_args()

    p_test = ParameterTesting(steps=args.steps)

    p_test.run_tests()
