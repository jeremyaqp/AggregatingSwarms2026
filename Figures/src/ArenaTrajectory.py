import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
from glob import glob

def plot_trajectory_with_velocity(x, y, cmap_name='viridis'):

    # Convert to numpy arrays
    x = np.asarray(x)
    y = np.asarray(y)
    
    # Sanity check
    if len(x) != len(y):
        raise ValueError("x and y must have the same length.")
    if len(x) < 2:
        raise ValueError("Trajectory must have at least two points.")

    # Compute instantaneous velocity (magnitude of dx/dt, dy/dt)
    # convolution of x and y to smooth the trajectory
    x = np.convolve(x, np.ones(30)/30, mode='valid')
    y = np.convolve(y, np.ones(30)/30, mode='valid')
    dx = np.diff(x)
    dy = np.diff(y)
    velocity = np.sqrt(dx**2 + dy**2) * 6.0  # assuming dt=1/6s
    velocity /= np.max(velocity)  # Normalize velocity for color mapping

    # Normalize velocity for colormap
    norm = plt.Normalize(vmin=velocity.min(), vmax=velocity.max())

    # Create colormap
    cmap = plt.get_cmap(cmap_name)
    print(cmap_name)
    fig, ax = plt.subplots()
    k=1
    # Plot segment by segment
    for i in range(0, len(velocity)-k, k):
        plt.plot(x[i:i+2+k], y[i:i+2+k], color=cmap(norm(velocity[i])), linewidth=2)
    plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, label='Normalized velocity')


if __name__ == "__main__":
    plt.rcParams.update({
        "text.usetex": True,
        "font.family": "Times New Roman",
        "font.size": 16
    })


    px_to_cm = 1.0 / 18.5
    x_padding = -270*px_to_cm

    files = glob("../../Data/ExampleArenaTrajectories.csv")
    for i, f in enumerate(files):
        title = f.split("/")[-1][:-9]
        df = pd.read_csv(f)
        bg = plt.imread("../../Data/ArenaTopView.png")
        for i in df["particle"].unique():
            if i >= 21 or i < 20:
                continue
            part = df[df["particle"] == i]
            part = part[part["frame"] < 6*1800]
            plot_trajectory_with_velocity(((part["x"]+20)*0.93)*px_to_cm, (2880-part["y"])*px_to_cm, cmap_name='seismic')
            print(i)
            if i%2 == 0:
                plt.imshow(bg, cmap="gray", extent=(x_padding, (3700 + x_padding)*px_to_cm, 0, (3000)*px_to_cm))
                plt.xlim(120*px_to_cm, (120+2830)*px_to_cm)
                plt.ylim(0*px_to_cm, (2830)*px_to_cm)
                plt.gca().set_aspect("equal")

                # center figure
                plt.xticks(ticks=np.linspace(120*px_to_cm, (120+2830)*px_to_cm, num=7), labels=[str(round(k)) for k in (np.linspace(0, 150.0, num=7)-75.0)])
                plt.yticks(ticks=np.linspace(0*px_to_cm, (2830)*px_to_cm, num=7), labels=[str(round(k)) for k in (np.linspace(0, 150.0, num=7)-75.0)])
                plt.xlabel(r"x (cm)")
                plt.ylabel(r"y (cm)")
                plt.savefig("../out/traj.png", dpi=300)
                plt.clf()
                exit()
            if i > 80:
                break