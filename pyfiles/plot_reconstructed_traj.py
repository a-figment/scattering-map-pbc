"""
@file plot_region_trajectories.py
@brief Reconstructs and plots billiard trajectories as determined from scattering map data.  
"""
import math
import numpy as np
import matplotlib
import matplotlib.animation as animation
import matplotlib.patches as patches
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib import colors

import trajectory_reconstructor as tr
import utils as myutils

def plot_channel(x,y,n_cells,ax):
    """@brief Plots a parallel polygonal billiard channel (`x`,`y`)
    @param n_cells Total number of cells either side of the E.C. at 0
    """
    n_cells = int(n_cells)
    ax.plot(x[0:3],y[0:3],color='black',linewidth=0.9)
    ax.plot(x[3:6],y[3:6],color='black',linewidth=0.9)
    if n_cells == 0:
        return ax
    if (n_cells % 2) != 0:
        n_cells -= 1
    if n_cells <= 1:
        n_cells = 2
    for i in range(1,int(n_cells/2)+1):
        ax.plot(x[0:3]+i,y[0:3],color='black',linewidth=0.9)
        ax.plot(x[3:6]+i,y[3:6],color='black',linewidth=0.9)
    for i in range(-int(n_cells/2),0):
        ax.plot(x[0:3]+i,y[0:3],color='black',linewidth=0.9)
        ax.plot(x[3:6]+i,y[3:6],color='black',linewidth=0.9)
    return ax

def plot_channel_filled(x,y,n_cells,ax):
    """@brief Plots channel by filling white space between boundaries
    @param n_cells Total number of cells either side of the E.C. at 0
    """
    if (n_cells % 2) == 0:
        n_cells += 1
    half_n = int(n_cells/2)
    x_lower = np.concatenate([np.tile(x,3) for x in range(-half_n, half_n+1)]) \
            + np.tile(x[0:3], n_cells)
    y_lower = np.tile(y[0:3], n_cells)    
    y_upper = np.tile(y[3:6], n_cells)    
    ax.fill_between(x_lower, y_lower, y_upper, color='gainsboro')
    ax.plot(x_lower, y_lower, color='black', linewidth=0.9)
    ax.plot(x_lower, y_upper, color='black', linewidth=0.9)

def plot_trajectory(Channel,Traj,n_cells):
    """@brief Plot a single trajectory 
    @param 
    """
    fig, ax = plt.subplots(dpi=300)
    ax.set_aspect(1)
    plot_channel_filled(Channel.x, Channel.y, n_cells,ax) 
    xt,yt = tr.reconstruct_lifted_trajectory(Channel,Traj)
    ax.set_xlim([-n_cells/2,n_cells/2])
    ax.plot(xt, yt, marker='o', markersize=0.05, linewidth=0.7, color='blue')  
    ax.scatter(xt[0],yt[0],color='green',s=1,marker='.')
    ax.scatter(xt[xt.size-1],yt[yt.size-1],color='green',s=5,marker='.')
    plt.show()

def plot_trajectories(Channel,Ensemble,N=5):
    """@brief Creates a 4x5 or 4x4 grid plot of individual trajectories
    - Reproduces Figs. 6 and 9

    TODO: Remove white space from figure and enable custom grid sizes
    TODO: Must exit/warn if data dimensions are exceed 4x5 or 4x4
    """
    if Ensemble.H.shape[0] < N:
        N = Ensemble.H.shape[0]

    # some random trajectories
    traj_indices = np.random.randint(0,Ensemble.H.shape[0],2*N)

    fig = plt.figure(figsize=(20,20),constrained_layout=True)
    gs = gridspec.GridSpec(2*N,2, figure=fig)
    for j in range(len(traj_indices)):
        s_traj = Ensemble.get_trajectory(traj_indices[j])
        xtraj,ytraj = tr.reconstruct_lifted_trajectory(Channel, s_traj)
        span = np.max(xtraj) - np.min(xtraj)
        if span > 10:
            print(f"Skipping trajectory with span {span}.")
        ax = plt.subplot(gs[j])
        ax.margins(x=0)
        ax.set_aspect(1)
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_xticks([])
        ax.set_yticks([])
        ax.plot(xtraj, ytraj, marker='o', markersize=0.05, linewidth=0.7, color='blue')  
        mpx = (np.max(xtraj)+np.min(xtraj)) / 2.
        plot_channel_filled(Channel.x+int(np.round(mpx)),Channel.y,int(np.round(span))+1,ax)
    gs.tight_layout(fig)
    plt.show()

def animate_tracer_trajectory(Channel,Ensemble):
    """@brief Tracer style animation frames of a single trajectory
    TODO: Everything...
          - Extend to many particles 
          - Plot with lines instead of points
          - Add colors, opacity of trailing etc.
    """
    x,y,l,e = tr.get_interpolated_cloud(Channel, Ensemble, s=0.05, min_frames=None)
    pid = np.random.randint(0,x.shape[0],1)
    lag = 30
    print(x.shape)
    def update(frame):
        start = max(0, frame - lag)
        #ax.set_xlim([x[frame,pid]-1, x[frame,pid]+1])
        pt.set_offsets(np.c_[x[start:frame,pid],y[start:frame,pid]])
        return pt,

    fig, ax = plt.subplots()
    ax.set_aspect(1)
    ax.set_xlim(-1, 1)
    plot_channel_filled(Channel.x,Channel.y,3,ax)
    pt = ax.scatter(x[0,pid], y[0,pid], marker='.', linewidths=0, s=36*(72./fig.dpi)**2)
    ani = animation.FuncAnimation(fig, update, frames=x.shape[0], interval=20, blit=True)
    plt.show()


def animate_particle_cloud(Channel,Ensemble,s=0.01,psize=None):
    """@brief Animates a particle cloud using matplotlib's `FuncAnimation`
    @param Channel 
    @param Ensemble 
    @param s Interpolation length 
    @param psize Particle size as psize * pixel
    """
    print("Reconstructing trajectories and interpolating ...")
    x_cloud,y_cloud,l_cloud = tr.get_interpolated_cloud(Channel,Ensemble,s)[0:3]
    print("Animating cloud...")

    # Animate
    fig, ax = plt.subplots(dpi=300)
    if psize is None:
        psize = 49*(72./fig.dpi)**2
    ax.set_facecolor('gray')
    ax.set_aspect(1)

    pt = ax.scatter(x_cloud[0], y_cloud[0], color=myutils.COLORS[l_cloud[0]], marker='.', linewidths=0, s=psize)
    pid = np.random.choice(range(Ensemble.H.shape[0]), 1)[0]  #index of particle to follow
    def update(frame):
        pt.set_offsets(np.c_[x_cloud[frame],y_cloud[frame]])
        pt.set_color(myutils.COLORS[l_cloud[frame]])
        #ax.set_xticks(range(int(x_cloud[pid,frame]-1), int(x_cloud[pid,frame])+2))  # Set x-ticks at integer values
        #ax.set_xlim([x_cloud[pid,frame]-1, x_cloud[pid,frame]+1])
        ax.set_xlim([-1.2,1.2])
        return pt,
    plot_channel_filled(Channel.x, Channel.y, 2*int(np.max(np.abs(x_cloud)))+2, ax)
    ani = animation.FuncAnimation(fig=fig, func=update, frames=x_cloud.shape[0], interval=30, blit=True)
    ani.save(myutils.RESULTS_PATH + 'cloud_animation.mp4', writer='ffmpeg')
    plt.show()
    print("Done!")


