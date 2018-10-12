import argparse
import datetime
import os
import subprocess

import numpy as np

from realistic_network import make_and_cd, BindingParameters


class LaunchQsub(object):
    def __init__(self, steps, lf=30, self_foreign=False, early_pos=False):
        self.steps = steps
        self.lf = lf

        self.self_foreign = self_foreign
        self.early_pos = early_pos
        self.simulation_name = "ODE_steps_" + str(self.steps)
        self.simulation_time = 10

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

        if self.self_foreign:
            if self.early_pos:
                q.write(
                    "python ~/SSC_python_modules/pysb_t_cell_network.py --steps {0} --lf {1} --early_pos_fb\n".format(
                        self.steps, self.lf))
            else:
                q.write("python ~/SSC_python_modules/pysb_t_cell_network.py --steps {0} --lf {1}\n".format(self.steps,
                                                                                                           self.lf))

        else:
            if self.early_pos:
                q.write(
                    "python ~/SSC_python_modules/pysb_t_cell_network.py --steps {0} --early_pos_fb".format(self.steps))
            else:
                q.write("python ~/SSC_python_modules/pysb_t_cell_network.py --steps {0}".format(self.steps))

        q.close()

    def launch(self):
        (stdout, stderr) = subprocess.Popen(["qsub {0}".format("qsub.sh")], shell=True, stdout=subprocess.PIPE,
                                            cwd=os.getcwd()).communicate()


class LaunchStochastic(LaunchQsub):

    def __init__(self, steps, lf=30, self_foreign=False):
        LaunchQsub.__init__(self, steps, lf=lf, self_foreign=self_foreign)

        self.simulation_name = "SSA_steps_" + str(self.steps)
        self.simulation_time = 20
        self.num_samples = 1000

        self.p_ligand = [int(i) for i in np.round(np.random.lognormal(6.0, 1.0, self.num_samples))]

    def make_qsub_script(self, ls, name):
        q = open("qsub.sh", "w")

        q.write("#PBS -m ae\n")
        q.write("#PBS -q short\n")
        q.write("#PBS -V\n")
        q.write("#PBS -l walltime={1},nodes=1:ppn=2 -N {0}\n\n".format(name, datetime.timedelta(
            minutes=self.simulation_time)))
        q.write("cd $PBS_O_WORKDIR\n")
        q.write("echo $PBS_JOBID > job_id\n\n")

        if self.self_foreign:
            q.write(
                "python ~/SSC_python_modules/pysb_t_cell_network.py --steps {0} --lf {1} --ssa {2}\n".format(self.steps,
                                                                                                             self.lf,
                                                                                                             ls))
        else:
            q.write("python ~/SSC_python_modules/pysb_t_cell_network.py --steps {0} --ssa {1}".format(self.steps, ls))

        q.close()

    def generate_qsub(self):
        home_directory = os.getcwd()
        sample = []
        for i in range(len(self.p_ligand)):
            directory = "sample_" + str(i)
            ls = self.p_ligand[i]

            sample.append(ls)
            simulation_name = self.simulation_name + "_" + str(i)

            os.makedirs(directory)
            print("Made " + directory)
            os.chdir(directory)
            print("Changed into directory: " + str(os.getcwd()))

            self.make_qsub_script(ls, simulation_name)
            (stdout, stderr) = subprocess.Popen(["qsub {0}".format("qsub.sh")], shell=True, stdout=subprocess.PIPE,
                                                cwd=os.getcwd()).communicate()
            os.chdir(home_directory)

        np.savetxt("Ligand_concentrations", sample, fmt='%f')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Submitting ode calculations as function of steps",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--steps', dest='steps', action='store', type=int, default=8,
                        help="number of KP steps.")
    parser.add_argument('--lf', dest='lf', action='store', type=int, default=30,
                        help="number of foreign ligands.")
    parser.add_argument('--ssa', dest='ssa', action='store_true', default=False,
                        help="flag for submitting stochastic calculations.")

    args = parser.parse_args()

    binding_parameters = BindingParameters()

    # if args.steps == 9:
    #     make_and_cd("{0}_step".format(args.steps))
    # else:
    make_and_cd("{0}_step".format(args.steps))

    sub_directories = ["Ls", "Ls_Lf_{0}".format(args.lf)]
    home_directory = os.getcwd()

    if args.ssa:
        sub_directories = ["Ls"]

    for directory in sub_directories:
        make_and_cd(directory)
        if directory == "Ls_Lf_{0}".format(args.lf):
            qsub = LaunchQsub(args.steps, lf=args.lf, self_foreign=True)
        else:
            qsub = LaunchQsub(args.steps)

        qsub.generate_qsub()
        qsub.launch()
        os.chdir(home_directory)
