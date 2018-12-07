import numpy as np


class PostProcessSingleMolecule(object):

    def __init__(self):
        self.column_names = open("column_names").readlines()[0].split()

        self.mean_trajectory = np.loadtxt("mean_traj")

        self.lf_indices = self.find_indices("Lf")
        self.ls_indices = self.find_indices("Ls")

        self.row_length = len(self.mean_trajectory[:, 0])

    def find_indices(self, string):
        indices = []
        for i in range(len(self.column_names)):
            if string in self.column_names[i]:
                indices.append(i)

        return indices

    def main(self):
        lf_output = np.zeros(self.row_length)
        for i in self.lf_indices:
            lf_output += self.mean_trajectory[:, i]

        ls_output = np.zeros(self.row_length)
        for i in self.ls_indices:
            ls_output += self.mean_trajectory[:, i]

        np.savetxt('total_output', np.c_[self.mean_trajectory[:, 0], self.mean_trajectory[:, 1],
                                         lf_output, ls_output], fmt='%1.3f')


if __name__ == "__main__":
    process = PostProcessSingleMolecule()

    process.main()
