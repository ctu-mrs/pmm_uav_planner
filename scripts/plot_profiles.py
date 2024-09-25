import sys
import string
import numpy as np
from csv import reader
import matplotlib.pyplot as plt
    
class pmm_sampled_trajectory_3d:
    def __init__(self) -> None:
        self.t = np.empty(0)
        self.p = np.empty((3,0))
        self.v = np.empty((3,0))
        self.a = np.empty((3,0))

    def load_trajectory_data_from_csv(self, file_name: str):
        with open(file_name, 'r') as file:
            csv_reader = reader(file)
            data = np.array([row for row in csv_reader], dtype=np.float64)
            
            data = data.T
            
            self.t = data[0, :]
            n_wp = self.t.shape[0]
            
            self.p = np.empty((3, n_wp))
            self.v = np.empty((3, n_wp))
            self.a = np.empty((3, n_wp))
            
            self.p = data[1:4, :]
            self.v = data[4:7, :]
            self.a = data[7:10, :]

    def plot_trajectory_data(self, color: list, axs, legend: list, title:string):
        if axs is None:
            # create new figure
            fig, axs = plt.subplots(3,1, figsize=(10,10))
            fig.suptitle(title, fontsize=16)

        marker_size = 2
        line_width = 2

        # position
        axs[0].plot(self.t, self.p[0,:],".-", color= color[0], linewidth=line_width, markersize=marker_size, label=legend[0])
        axs[0].plot(self.t, self.p[1,:],".-", color= color[1], linewidth=line_width, markersize=marker_size, label=legend[1])
        axs[0].plot(self.t, self.p[2,:],".-", color= color[2], linewidth=line_width, markersize=marker_size, label=legend[2])
        axs[0].set_xlabel("t [s]")
        axs[0].set_ylabel("p [m]")
        axs[0].grid(True)
        axs[0].legend()

        # velocity
        axs[1].plot(self.t, self.v[0,:], ".-", color= color[0], linewidth=line_width, markersize=marker_size)
        axs[1].plot(self.t, self.v[1,:], ".-", color= color[1],linewidth=line_width, markersize=marker_size)
        axs[1].plot(self.t, self.v[2,:], ".-", color= color[2], linewidth=line_width, markersize=marker_size)
        axs[1].set_xlabel("t [s]")
        axs[1].set_ylabel("v [m/s]")
        axs[1].grid(True)

        # acceleration
        axs[2].plot(self.t, self.a[0,:], ".-", color= color[0], linewidth=line_width, markersize=marker_size)
        axs[2].plot(self.t, self.a[1,:], ".-", color= color[1],linewidth=line_width, markersize=marker_size)
        axs[2].plot(self.t, self.a[2,:], ".-", color= color[2], linewidth=line_width, markersize=marker_size)
        axs[2].set_xlabel("t [s]")
        axs[2].set_ylabel("a [m/s^2]")
        axs[2].grid(True)
        return axs

def visualize_sampled_pmm_mp_3d(trajectory_file):

    tr_optim2 = pmm_sampled_trajectory_3d()
    tr_optim2.load_trajectory_data_from_csv(trajectory_file)

    axs = tr_optim2.plot_trajectory_data(['r', 'g', 'b'], None, ['x', 'y', 'z'], 'Trajectory Profiles')

if __name__=="__main__":

    if len(sys.argv) < 2:
        print("ERROR: Please specify trajectory file!")

    print("Visualizing Trajectory Profile")
    visualize_sampled_pmm_mp_3d(sys.argv[1])

    plt.show()