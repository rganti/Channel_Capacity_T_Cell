import matplotlib.pyplot as plt
import numpy as np

from post_process import load, getcolumnames


class PlotTrajectory(object):
    def __init__(self, directory="./"):
        self.directory = directory
        self.column_names = getcolumnames(load(self.directory + "column_names"))

        self.mean_trajectory = np.loadtxt(self.directory + "mean_traj")
        self.time = self.mean_trajectory[:, 0]

    def plot_trajectory(self):
        plt.plot(self.time, self.mean_trajectory[:, -1], label=self.column_names[-1])
        plt.plot(self.time, self.mean_trajectory[:, -2], label=self.column_names[-2])

        plt.legend()
