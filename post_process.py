#!/usr/bin/python
import argparse
import os
import re

import numpy as np


def load(filename):
    # print("Loading data...")
    data = list(map(lambda x: x, open(str(filename)).readlines()))
    return data


def getcolumnames(data):
    columnames = re.sub('[#\[/\]]', '', data[0]).split()
    return columnames


class PostProcess(object):
    def __init__(self, trajectory_file):
        self.trajectory_file = trajectory_file
        self.data = load(self.trajectory_file)
        self.column_names = getcolumnames(self.data)

    def write_columns(self):
        f = open("column_names", "w")
        for item in self.column_names:
            f.write("{0} ".format(item))
        f.write("\n")
        f.close()

    def insert_hash(self):
        self.data[0] = "# " + self.data[0]

    def write_data(self):
        f = open("hashed_{0}".format(self.trajectory_file), "w")
        f.write("".join(map(lambda x: x, self.data)))
        f.close()
        # os.remove(self.trajectory_file)

    def main(self):
        self.insert_hash()
        self.write_data()


def average_trajectories(num_files, run_time, time_step):
    all_data = []
    num_files = len([name for name in os.listdir(".") if "hashed" in name])
    print(num_files)
    for i in range(1, num_files):
        array = np.loadtxt("hashed_traj_{0}".format(i))
        print(array.shape)
        if float(array.shape[0]) == run_time/time_step:
            all_data += [array]
        elif run_time/time_step == 1:
            all_data += [array]
        os.remove("hashed_traj_{0}".format(i))

    np_all_data = np.array(all_data)
    print(np_all_data.shape)
    # all_data = np.array([np.loadtxt("hashed_traj_{0}".format(i)) for i in range(1, num_files)])
    # print(all_data.shape)
    average_trajectory = np.mean(all_data, axis=0)
    np.savetxt("mean_traj", average_trajectory, fmt="%f")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Args for post-processing KP output.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--num_files', dest='num_files', action='store', type=int,
                        help="number of trajectory files generated.")
    parser.add_argument('--run_time', dest='run_time', action='store', type=int,
                        help="run time for trajectories.")
    parser.add_argument('--time_step', dest='time_step', action='store', type=float,
                        help="time step used for simulations.")

    args = parser.parse_args()

    # num_files = args.num_files
    run_time = args.run_time
    time_step = args.time_step

    post_process = PostProcess("traj_1")
    post_process.write_columns()

    num_files = len([name for name in os.listdir(".") if "traj" in name])
    for i in range(num_files - 1):
       post_process = PostProcess("traj_{0}".format(i+1))
       post_process.main()

    average_trajectories(num_files, run_time, time_step)
