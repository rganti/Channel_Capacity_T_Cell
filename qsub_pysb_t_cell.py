import argparse
import datetime
import os
import subprocess

from realistic_network import make_and_cd, BindingParameters


class LaunchQsub(object):
    def __init__(self, steps, lf=30, self_foreign=False):
        self.steps = steps
        self.lf = lf

        self.self_foreign = self_foreign
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


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Submitting ode calculations as function of steps",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--steps', dest='steps', action='store', type=int, default=8,
                        help="number of KP steps.")
    parser.add_argument('--lf', dest='lf', action='store', type=int, default=30,
                        help="number of foreign ligands.")

    args = parser.parse_args()

    binding_parameters = BindingParameters()

    if args.steps == 9:
        make_and_cd("{0}_step_k_pos_{1}".format(args.steps, binding_parameters.k_positive_loop))
        # make_and_cd("{0}_step".format(args.steps))
    else:
        make_and_cd("{0}_step".format(args.steps))

    sub_directories = ["Ls", "Ls_Lf_{0}".format(args.lf)]
    home_directory = os.getcwd()

    for directory in sub_directories:
        make_and_cd(directory)
        if directory == "Ls_Lf_{0}".format(args.lf):
            qsub = LaunchQsub(args.steps, lf=args.lf, self_foreign=True)
        else:
            qsub = LaunchQsub(args.steps)

        qsub.generate_qsub()
        qsub.launch()
        os.chdir(home_directory)
