import sys
import numpy as np
import matplotlib.pyplot as plt
import csv

def load_trajectory_samples_pmm_row(file, header=False):
    states = []
    with open(file, 'r') as csvfile:
        csvreader = csv.reader(csvfile)
        if header:
            header = next(csvreader)
        for row in csvreader:
            col = []
            for c in row:
                col.append(float(c))
            states.append(col)
    return np.array(states)

def plot_3d_positions_graph(pmm_samples):
    fig2 = plt.figure(figsize=(12, 7))
    ax3d = fig2.add_subplot(111, projection='3d', computed_zorder=False)

    velocity_norms = np.sqrt(pmm_samples[:, 4] ** 2 + pmm_samples[:, 5] ** 2 + pmm_samples[:, 6] ** 2)
    velocities_plot = ax3d.scatter(pmm_samples[:, 1], pmm_samples[:, 2], pmm_samples[:, 3], c=velocity_norms, cmap='jet', s=0.35, zorder=-1)

    ax3d.axis('equal')
    ax3d.set_xlabel('x [m]', labelpad=10)
    ax3d.set_ylabel('y [m]', labelpad=10)
    ax3d.set_zlabel('z [m]', labelpad=10)
    ax3d.grid(False)
    fig2.colorbar(velocities_plot, label="velocity [m/s]", pad=0.15,
                  ticks=[min(velocity_norms) + (max(velocity_norms) - min(velocity_norms)) * i / 10 for i in range(0, 11, 2)], shrink=0.7)
    fig2.subplots_adjust(left=0.05, right=0.95, bottom=0.05, top=0.95, wspace=0.05, hspace=0.2)


if __name__=="__main__":
    
    if len(sys.argv) < 2:
        print("ERROR: Please specify trajectory file!")

    print("Visualizing Trajectory")

    tr_samples = load_trajectory_samples_pmm_row(sys.argv[1])
    plot_3d_positions_graph(tr_samples)

    plt.show()