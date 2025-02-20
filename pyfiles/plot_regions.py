"""
@brief Plots the coordinate space using rejection sampled numerical data. 
"""

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

import utils as myutils

def animate_regions(Channel, H, Theta):
    """@brief Animates the `evolution' of the coordinate space 
    @param H 
    @param Theta
    """


def plot_coordinate_space_region(d, rlabels, H_FILE, THETA_FILE):
    """@brief Plots regions of the coordinate space
    `H` and `Theta` are 2D arrays of 
    @param d Channel width
    @param H_FILE Filename for heights 
    @param THETA_FILE '' angles
    """
    xlimits = [0.,d]
    aspect = 0.25
    if d > .9:
        aspect = 0.5

    fig, ax = plt.subplots(dpi=150)
    ax.set_aspect(aspect) 
    ax.set_xlabel('h')
    ax.set_ylabel('\u03B8')
    ax.set_ylim([-np.pi/2, np.pi/2])
    ax.set_xlim(xlimits)
    with open(H_FILE, "r") as h, open(THETA_FILE, "r") as t:   # Arrays are ragged if using weights
        for i, (hline, tline) in enumerate(zip(h,t)):
            H     = np.fromiter(map(np.float64, hline.split()), dtype=np.float64)
            Theta = np.fromiter(map(np.float64, tline.split()), dtype=np.float64)
            ax.plot(H, Theta, color=myutils.COLORS[i], marker=',', markersize=0.2, linestyle='')
            ax.text(x=np.mean(H), y=np.mean(Theta), s=str(rlabels[i]), fontsize=5)

    ax.set_title(f"Coordinate space (d = {d})")
    plt.show()
    fig.savefig("../results/d-" + str(d) + "_coordinate_space.png", dpi=300)
    ax.cla()
    ax.clear()
    plt.close(fig)

def save_coordinate_space_plots():
    """@brief Saves coordinate space plots
    """
    # Plot static regions for each study 
    studies = [myutils.STUDY_PREFIX + str(i) for i in [0,1]]
    for study in studies:
        params = myutils.get_study_parameters(study)
        file_dict = myutils.get_study_files(study, warn=True)
        if Path(file_dict["Regions-H"]).exists() and Path(file_dict["Regions-Theta"]).exists():
            plot_coordinate_space_region(params["width"], params["region_labels"], file_dict["Regions-H"], file_dict["Regions-Theta"])

if __name__ == '__main__':
    save_coordinate_space_plots()



