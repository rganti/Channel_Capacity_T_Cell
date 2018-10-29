import argparse
import datetime
import os
import pickle
import subprocess

import numpy as np
import pandas as pd

from realistic_network import make_and_cd, BindingParameters


class ParameterTesting(object):
    def __init__(self, steps=8):
        self.binding_constants = BindingParameters()

        self.lf = 30
        self.sub_directories = ["Ls", "Ls_Lf_{0}".format(self.lf)]
        self.steps = steps
        #
        # self.k_8_1 = np.round(np.logspace(-3.0, -2.5, num=10), 6)
        # self.param_grid = {'k_8_1': self.k_8_1}
        #
        # if self.steps > 8:
        #     self.param_grid['k_9_1'] = np.round(np.logspace(-5.0, -2.5, num=8), 6)  # np.linspace(0.00005, 0.0003, 3)
        #
        # self.grid = ParameterGrid(self.param_grid)
        # print(list(self.grid))

        self.paths = []
        self.parameters = []
        self.home_directory = os.getcwd()

        self.simulation_name = "ODE_steps_" + str(self.steps)
        self.simulation_time = 10

    def generate_qsub(self, self_foreign=False):

        q = open("qsub.sh", "w")
        q.write("#PBS -m ae\n")
        q.write("#PBS -q short\n")
        q.write("#PBS -V\n")
        q.write("#PBS -l walltime={1},nodes=1:ppn=2 -N {0}\n\n".format(self.simulation_name,
                                                                       datetime.timedelta(
                                                                           minutes=self.simulation_time)))
        q.write("cd $PBS_O_WORKDIR\n")
        q.write("echo $PBS_JOBID > job_id\n\n")

        if self_foreign:
            q.write(
                "python ~/SSC_python_modules/pysb_t_cell_network.py --steps {0} --lf {1}\n".format(
                    self.steps,
                    self.lf))
        else:
            q.write("python ~/SSC_python_modules/pysb_t_cell_network.py --steps {0}".format(self.steps))
        q.close()

    def launch(self):
        (stdout, stderr) = subprocess.Popen(["qsub {0}".format("qsub.sh")], shell=True, stdout=subprocess.PIPE,
                                            cwd=os.getcwd()).communicate()

    def make_launch_simulations(self, set):
        home = os.getcwd()
        for directory in self.sub_directories:
            make_and_cd(directory)
            pickle_out = open("parameters.pickle", "wb")
            pickle.dump(set, pickle_out)
            pickle_out.close()

            if directory == "Ls_Lf_{0}".format(self.lf):
                self.generate_qsub(self_foreign=True)
            else:
                self.generate_qsub()

            if args.run:
                self.launch()
            os.chdir(home)

    def create_submit(self, count, param_grid):
        file_path = "{0}_step_{1}".format(self.steps, count)
        self.paths.append(file_path)
        print("param_grid " + str(param_grid))
        self.parameters.append(str(param_grid))

        make_and_cd(file_path)
        self.make_launch_simulations(param_grid)
        os.chdir(self.home_directory)

    def run_parameter_search(self):
        count = 0

        on_rates = np.linspace(0.5, 5.0, 10)
        discount = [1.0, 5.0, 10.0, 50.0, 100.0]

        param_grid = {}

        for on_rate in on_rates:
            param_grid['k_lck_on_RL'] = on_rate / self.binding_constants.initial.lck_0
            param_grid['k_p_on_R_pmhc'] = on_rate
            param_grid['k_zap_on_R_pmhc'] = on_rate / self.binding_constants.initial.zap_0
            param_grid['k_p_on_zap_species'] = on_rate
            param_grid['k_lat_on_species'] = on_rate / self.binding_constants.initial.lat_0
            param_grid['kp_on_lat1'] = on_rate

            if self.steps > 7:
                for i in discount:
                    param_grid['kp_on_lat2'] = on_rate / i
                    self.create_submit(count, param_grid)

                    count += 1
            else:
                self.create_submit(count, param_grid)
                count += 1

        df = pd.DataFrame({'file_path': self.paths})
        df.to_csv("./file_paths", sep='\t')

        df_2 = pd.DataFrame({'file_path': self.paths, 'parameters': self.parameters})
        df_2.to_csv("./parameters", sep='\t')

        pickle_out = open("parameter_range.pickle", "wb")
        pickle.dump(param_grid, pickle_out)
        pickle_out.close()

    # def run_tests(self):
    #     paths = []
    #     parameters = []
    #
    #     count = 0
    #     home_directory = os.getcwd()
    #
    #     for set in list(self.grid):
    #         file_path = "{0}_step_{1}".format(self.steps, count)
    #         paths.append(file_path)
    #         parameters.append(str(set))
    #
    #         make_and_cd(file_path)
    #         self.make_launch_simulations(set)
    #         os.chdir(home_directory)
    #
    #         count += 1
    #
    #     df = pd.DataFrame({'file_path': paths})
    #     df.to_csv("./file_paths", sep='\t')
    #
    #     df_2 = pd.DataFrame({'file_path': paths, 'parameters': parameters})
    #     df_2.to_csv("./parameters", sep='\t')
    #
    #     pickle_out = open("parameter_range.pickle", "wb")
    #     pickle.dump(self.param_grid, pickle_out)
    #     pickle_out.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Submitting ode calculations as function of steps",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--steps', dest='steps', action='store', type=int, default=9, help="number of KP steps.")
    parser.add_argument('--run', action='store_true', default=False, help='Flag for submitting simulations.')

    args = parser.parse_args()

    p_test = ParameterTesting(steps=args.steps)

    p_test.run_parameter_search()
