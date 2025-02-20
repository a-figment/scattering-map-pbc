"""
@brief Free functions used for reconstructing trajectories from scattering map data

NOTE 1: Itineraries from Tables 1,2,3 must be adjusted to account for the additional
        edges. e.g. \\mathcal{L}\\Gamma_2\\Gamma_1\\mathcal{R} = 5 31 2

NOTE 2: The routines are intended for visualisation. If needed, they could be optimised 
        by calculating all tilings and intersections in advance
"""

import math
import numpy as np
#import utils as myutils

def is_close(a, b, rel_tol=1e-15, abs_tol=0.0):
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

def intersect_ray_with_edge(x, y, h, theta, edge_id):
    """@brief Finds the intersection between ray (`h`,`theta`) and line segment

    NOTE: Always exists for valid itineraries/unfoldings
    @param edge_id of polygon (`x`,`y`)
    @return x_soln,y_soln Intersection point on the unfolded surf
    """
    m = (y[edge_id+1]-y[edge_id]) / (x[edge_id+1] - x[edge_id])
    x_soln = (h - y[edge_id] + m*x[edge_id]) / (m - math.tan(theta))
    y_soln = h + x_soln * math.tan(theta) 
    return x_soln, y_soln

def reflect_about_edge(x, y, edge_id):
    """@brief Reflects a polygon (`x`,`y`) about its edge `edge_id`. 
    @return x_Refl,y_Refl vertices of the reflected polygon
    """
    x_refl = np.copy(x)
    y_refl = np.copy(y)

    # Reflections about a vertical edge
    if (is_close(x[edge_id],x[edge_id+1])):
        for i in range(len(x)):
            if (i == edge_id or i == edge_id+1):
                continue
            dx = x[edge_id] - x[i]
            x_refl[i] = x[i] + 2.*dx
    else:
        edge_slope = (y[edge_id+1] - y[edge_id])/(x[edge_id+1] - x[edge_id])
        edge_intcpt = (x[edge_id+1] * y[edge_id] - x[edge_id] * y[edge_id+1]) / (x[edge_id+1] - x[edge_id])
        for i in range(len(x)):
            if (i == edge_id or i == edge_id+1):
                continue
            d = (x[i] + (y[i] - edge_intcpt)*edge_slope)/(1. + edge_slope**2)
            x_refl[i] = 2.*d - x[i]
            y_refl[i] = 2.*d*edge_slope - y[i] + 2.*edge_intcpt
    return x_refl, y_refl

def get_point_from_arc_from_point(Channel, x_refl, y_refl, x_itsct, y_itsct, edge_id):
    """@brief Gets collision point on the elementary cell from the point on the unfolded channel

    NOTE: Edge labels are preserved on reflection, despite the alternation between 
    clockwise and counterclockwise vertex ordering (which is not relevant here). 

    @param Channel The elementary cell
    @param x_refl, y_refl Coordinates of a reflected `Channel`
    @param x_itsct, y_itsct A point on the reflected `Channel`'s boundary
    @param edge_id The segment on which the point lies
    @return Collision point on the E.C
    """
    # Get arc length from intersection point
    dx, dy = x_refl[edge_id+1] - x_refl[edge_id], y_refl[edge_id+1] - y_refl[edge_id]
    t = ((x_itsct - x_refl[edge_id]) * dx + (y_itsct - y_refl[edge_id]) * dy) / (dx**2 + dy**2)
    arc_len = Channel.cum_len[edge_id] + t * Channel.edge_len[edge_id]

    # Get intersection point from arc length on the base cell
    u = (arc_len - Channel.cum_len[edge_id]) / Channel.edge_len[edge_id]
    x_result = Channel.x[edge_id] + u * (Channel.x[edge_id+1] - Channel.x[edge_id])
    y_result = Channel.y[edge_id] + u * (Channel.y[edge_id+1] - Channel.y[edge_id])

    return x_result, y_result

def reconstruct_trajectory(Channel, h, theta, itinerary):
    """@brief Reconstructs a trajectory within the elementary cell according to `itinerary`. 

    For each character in `itinerary`, an unfolding of `Channel` is produced together with 
    the point of intersection between the ray (`h`,`theta`) and `Channel`, (`x_itsct`, `y_itsct`). 
    This point is then mapped back to the elementary cell, completing the reconstruction

    @param h Singleton scattering coordinate
    @param theta Singleton scattering coordinate in (-pi/2, pi/2)
    @param Channel The elementary cell
    @param itinerary Ordered list of integers
    @return x_soln,y_soln Reconstructed trajectory coordinates
    """
    x_soln, y_soln = np.empty(len(itinerary)), np.empty(len(itinerary))
    x_refl, y_refl = np.copy(Channel.x), np.copy(Channel.y) 
    for i in range(len(itinerary)):
        edge_id = itinerary[i]
        # NOTE: Iterative procedure - previous x_refl,y_refl must be supplied 
        x_refl, y_refl = reflect_about_edge(x_refl, y_refl, edge_id) 
        x_itsct, y_itsct = intersect_ray_with_edge(x_refl, y_refl, h, theta, edge_id) 
        x_soln[i], y_soln[i] = get_point_from_arc_from_point(Channel, x_refl, y_refl, x_itsct, y_itsct, edge_id)
    return x_soln, y_soln

def reconstruct_lifted_trajectory(Channel, Traj):
    """@brief Converts itinerary list to lifted billiard trajectory 

    @param Channel The elementary cell
    @param Traj An ensemble dataclass of 1D arrays characterising scattering map data
    @return xtraj,ytraj coordinate arrays corresponding to trajectory path
    TODO: Add an 'up_to_time' argument for reconstructions up to a traj length?
    """

    # Number of collisions ~ number of chars in itinerary
    total_chars = sum(len(s) for s in Traj.Itineraries)
    num_collisions = total_chars - len(Traj.Itineraries) - len(Traj.Itineraries[-1][1:-1])

    # Prep angles - reconstruction performed for ancoming angles (-pi/2, pi/2)
    thetas = np.copy(Traj.Theta)
    lhs_idx = np.where(~((thetas > -np.pi/2.) & (thetas < np.pi/2.)))[0] 
    thetas[lhs_idx] = np.pi - thetas[lhs_idx]
    symmetry_coeff = np.ones(len(thetas), dtype=int)
    symmetry_coeff[lhs_idx] = -1*symmetry_coeff[lhs_idx] # used to mirror x-coords

    xtraj, ytraj = np.empty(num_collisions), np.empty(num_collisions)
    xtraj[0],ytraj[0] = Traj.Positions[0], Traj.H[0]
    count = 1
    for i in range(len(Traj.Itineraries)-1):
        itin_chars = list(Traj.Itineraries[i]) # e.g. ['L','0','3', 'L']
        itinerary = list(map(int,itin_chars[1:(len(itin_chars)-1)])) # [0, 3]

        # Trajectory coordinates on the E.C.
        #xt, yt = np.empty(len(itinerary)), np.empty(len(itinerary))
        #if (Traj.Theta[i] > -np.pi/2. and Traj.Theta[i] < np.pi/2.):
        #    xt, yt = reconstruct_trajectory(Channel,Traj.H[i],Traj.Theta[i],itinerary)
        #else:
        #    xt, yt = reconstruct_trajectory(Channel,Traj.H[i],np.pi-Traj.Theta[i],itinerary)
        #    xt = -xt # mirror trajectory
        
        # TODO: 
        xt, yt = reconstruct_trajectory(Channel,Traj.H[i],thetas[i],itinerary)
        xt = symmetry_coeff[i] * xt

        xtraj[range(count, count+len(xt)+1)] = np.insert(xt+Traj.Positions[i], [xt.size], [Traj.Positions[i+1]])
        ytraj[range(count, count+len(xt)+1)] = np.insert(yt,                   [yt.size], [Traj.H[i+1]])
        count += (len(xt) + 1)
    return xtraj, ytraj

#def get_labels_from_interpolated_trajectory(x_intp, labels):
#    """@brief Maps the interpolated particle position to its label
#    
#    Example: [0,1/4,3/4,1]->[6,6,6,6]
#    """
#    x_trunc = np.trunc(x_intp) 
#    max_x = np.max(x_trunc)
#    # Cells [-1,0],[0,1] detected with np.signbit. 
#    signed_zeros = (x_trunc.astype(int) == 0) & np.signbit(x_trunc) 
#    unsigned_zeros = (x_trunc.astype(int) == 0) & (np.signbit(x_trunc) == False) 
#    # Only need to edge detect with x_trunc - replace with any unique element 
#    x_trunc[signed_zeros] = max_x + 200
#    x_trunc[unsigned_zeros] = max_x + 300
#    all_crossover = (x_trunc[:-1] != x_trunc[1:])  # Crossing events
#    crossing_locs = np.where(all_crossover)[0] + 1 # Crossing indices
#    crossing_bits = np.zeros(len(x_trunc), dtype=int)
#    # split and replicate labels between each crossing event
#    split_indices = np.split(crossing_bits, crossing_locs)
#    result = np.concatenate([np.full_like(segment, labels[i]) for i, segment in enumerate(split_indices[0:(len(labels))])])
#    return result

#def get_interpolated_trajectory(xr_traj, yr_traj, labels, thetas, heights, mf, s):
def get_interpolated_trajectory(xr_traj, yr_traj, labels, mf, s):
    """@brief Divides reconstructed trajectory into segments of length `s` (accounting for endpoints). 

    @return ends Sequence of indices that determine the range of frames for which
    a particle is within a given cell 

    NOTE: Dwell times have no lower bound. If interpolated points are at a fixed
    distance apart from eachother along the polyline, it is possible to not detect
    a crossing event, meaning that a change in the label is missed. We force an update
    on the interpolated label simply by checking if the next point on the reconstructed
    trajectory lies on \\mathcal{L} or \\mathcal{R}

    TODO: too slow, unravel loop
    TODO: only the 'ends' are needed to reconstruct     
    """
    fp = np.sqrt(np.diff(xr_traj)**2 + np.diff(yr_traj)**2)
    total_len = np.sum(fp)

    # directions long each seg
    x_dir = (xr_traj[1:] - xr_traj[:-1])/fp
    y_dir = (yr_traj[1:] - yr_traj[:-1])/fp
    total_elem = int(total_len // s) + 1

    xt, yt = np.empty(total_elem, dtype=np.float64), np.empty(total_elem, dtype=np.float64)
    lt = np.empty(total_elem, dtype=int) # labels
    ends = -2*np.ones(total_elem, dtype=int) # padded 'ends'
    xt[0],yt[0],lt[0],ends[0] = xr_traj[0], yr_traj[0], labels[0], 0

    lcount = 0
    leftover = 0.
    start = 1
    for i in range(len(xr_traj)-1):
        num_segs = int((fp[i]+leftover)//s)+1 # number of segments that fit into current seg
        I = np.arange(1, num_segs)
        # Fill
        end = start+len(I)
        xt[start:end] = xr_traj[i] + (I*s-leftover) * x_dir[i]
        yt[start:end] = yr_traj[i] + (I*s-leftover) * y_dir[i]
        lt[start:end] = np.repeat(labels[lcount], end-start)
        start = end
        
        # Next point on reconstructed polyline is on an entrance/exit line
        # TODO: np.where(xr_traj[1:] == xr_traj[1:].astype(int))...
        if xr_traj[i+1] == xr_traj[i+1].astype(int):
            lcount += 1
            ends[lcount] = end

        if num_segs == 1:       # Interpolation step exceeds segment length
            leftover += fp[i]
        else:                   # Get dist between reconstructed point and last interp. point
            leftover = np.sqrt((xt[end-1] - xr_traj[i+1])**2 + (yt[end-1] - yr_traj[i+1])**2)

        if end > mf+1:
            return xt,yt,lt,ends
    return xt,yt,lt,ends

def get_interpolated_cloud(Channel,Ensemble,s=0.01, min_frames=None):
    """@brief Prepares an `Ensemble` for animation by interpolating along each trajectory

    @param Ensemble dataclass of NDArray's used for holding scattering map data
    @param s Interpolation length
    @return Reconstructed and interpolated trajectory data

    TODO: The size of `Ensemble.get_trajectory(i)` should reflect `min_frames` (redundant calcs)
    """

    # No. frames determined out of scope for incremental processing
    if min_frames is None:
        # No. collisions in each reconstructed trajectory
        r_traj_length = np.vectorize(len)(Ensemble.Itineraries).sum(axis=1) # Total number of chars
        r_traj_length = r_traj_length - Ensemble.H.shape[1] - [len(x[1:-1]) for x in Ensemble.Itineraries[:,-1]]

        # Shortest trajectory determines minimum number of frames to animate for
        min_traj_length = np.min(np.sum(Ensemble.Tau, axis=1)) # shortest trajectory length
        min_frames = int(min_traj_length//s) #+ np.min(r_traj_length)

    interp_traj_x = np.empty((min_frames,Ensemble.H.shape[0]), dtype=np.float64)
    interp_traj_y = np.empty((min_frames,Ensemble.H.shape[0]), dtype=np.float64)
    interp_traj_l = np.empty((min_frames,Ensemble.H.shape[0]), dtype=int) 
    interp_traj_e = -2*np.ones((min_frames,Ensemble.H.shape[0]), dtype=int) 
    K = slice(0, min_frames)
    # Reconstruct and interpolate trajectories 
    for i in range(Ensemble.H.shape[0]): 
        Traj = Ensemble.get_trajectory(i)
        xr_traj,yr_traj = reconstruct_lifted_trajectory(Channel, Traj)
        X,Y,L,E = get_interpolated_trajectory(xr_traj, yr_traj, Traj.Labels, min_frames, s)
        interp_traj_x[:,i],interp_traj_y[:,i],interp_traj_l[:,i],interp_traj_e[:,i] = X[K], Y[K], L[K], E[K]
    return interp_traj_x, interp_traj_y, interp_traj_l, interp_traj_e



