import matplotlib.pyplot as plt

import utils as myutils                                # <-- Wrapper functions for reading config.json / loading data
import trajectory_reconstructor as tr                  # <-- Converts scattering map trajectory into billiard trajectory
import plot_reconstructed_traj as prt                  # <-- Misc. trajectory plotting / animating routines
import plot_regions as pr                              # <-- Plots and animates the coordinate space



study    = myutils.STUDY_PREFIX + "0"                  # <-- Specify study/experiment 
Channel  = myutils.get_polygonal_channel(study)        # <-- Load channel dataclass 
Ensemble = myutils.get_ensemble(study)                 # <-- Load ensemble dataclass 

xp = myutils.get_trajectory_member(study, "Positions") # <-- Equiv. to `Ensemble.Positions`
print(Ensemble)
Traj = Ensemble.get_trajectory(2)                      # <-- Get a trajectory from the ensemble
prt.plot_trajectory(Channel,Traj,n_cells=5)            # <-- Plot trajectory no. 2 (`Traj`)
prt.plot_trajectories(Channel,Ensemble,N=3)            # <-- Plot 2*N random trajectories
pr.save_coordinate_space_plots()                       # <-- Saves/shows image of cspace points
prt.animate_tracer_trajectory(Channel, Ensemble)       # <-- Tracer style animation of a trajectory 
prt.animate_particle_cloud(Channel,Ensemble,s=0.01,psize=None)

