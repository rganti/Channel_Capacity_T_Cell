import argparse
import os
import subprocess

from post_process import load


def restart_simulations(run=False):
    home_directory = os.getcwd()
    count = 0
    for i in range(1000):
        os.chdir("sample_{0}/".format(i))
        if not os.path.exists("mean_traj"):
            print(str(os.getcwd()))
            count += 1

            if run:
                qsub = load("qsub.sh")
                qsub[3] = qsub[3].replace("0:15:00", "0:45:00")
                qsub_file = open("qsub.sh", "w")
                for item in qsub:
                    qsub_file.write("{0}".format(item))
                qsub_file.close()

                (stdout, stderr) = subprocess.Popen(["qsub {0}".format("qsub.sh")], shell=True, stdout=subprocess.PIPE,
                                                    cwd=os.getcwd()).communicate()
        os.chdir(home_directory)

    print(str(count))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Submitting job for calculating P(C0) as function of steps",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--run', action='store_true', default=False,
                        help='Flag for submitting simulations.')

    args = parser.parse_args()

    restart_simulations(run=args.run)
