import datetime
import os
import subprocess


def generate_qsub(simulation_name, simulation_time):
    q = open("qsub.sh", "w")
    q.write("#PBS -m ae\n")
    q.write("#PBS -q short\n")
    q.write("#PBS -V\n")
    q.write("#PBS -l walltime={1},nodes=1:ppn=1 -N {0}\n\n".format(simulation_name,
                                                                   datetime.timedelta(minutes=simulation_time)))
    q.write("cd $PBS_O_WORKDIR\n\n")
    q.write("echo $PBS_JOBID > job_id\n")
    q.write("python ~/SSC_python_modules/pysb_t_cell_network.py \n")
    q.close()


def main_script():
    (stdout, stderr) = subprocess.Popen(["qsub {0}".format("qsub.sh")], shell=True, stdout=subprocess.PIPE,
                                        cwd=os.getcwd()).communicate()


generate_qsub("Ode_negative_fb", 1)
main_script()
