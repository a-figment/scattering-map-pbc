import math
import json
from pathlib import Path
import numpy as np
from dataclasses import dataclass, field, fields
from numpy.typing import NDArray
import pandas as pd

# Project paths and dirs (TODO: set defaults in config.json)
LOC = str(Path(__file__).resolve().parent) # TODO: 
FILE_JSON    = "../config.json"
STUDY_PREFIX = "study_"
DATA_PATH    = "../data/" + STUDY_PREFIX
RESULTS_PATH = "../results/"

COLORS = np.array([
    (0., 128/255., 1.),             # {0,128,255.}
    (0., 166/255., 90/255.),        # {0,166,90}
    (1., 153/255., 0.),             # {255.,153,0}
    (102/255., 51/255., 153/255.),  # {102,51,153}
    (0., 166/255., 90/255.),        # {0,166,90}
    (232/255., 69/255., 115/255.),  # {232,69,115}
    (162/255., 116/255., 50/255.),  # {162,116,50}
    (0., 128/255., 128/255.),       # {0,128,128}
    (128/255., 128/255., 128/255.), # {128,128,128}
    (1., 215/255., 0.),             # {255.,215,0}
    (31/255., 119/255., 180/255.),  # {31,119,180}
    (214/255., 39/255., 40/255.),   # {214,39,40}
    (1., 187/255., 120/255.),       # {255.,187,120}
    (148/255., 103/255., 189/255.), # {148,103,189}
    (222/255., 111/255., 161/255.), # {222,111,161}
    (204/255., 184/255., 114/255.), # {204,184,114}
    (130/255., 190/255., 180/255.), # {130,190,180}
    (191/255., 191/255., 191/255.), # {191,191,191}
    (1., 152/255., 150/255.),       # {1,152,150}
    (1., 242/255., 0.)              # {1,242,0}
])

@dataclass
class Ensemble:
    """@brief Collection of arrays that hold scattering map data
    All members are NxM with N = no. particles and M = no. iterations or time steps
    """
    H: NDArray[np.float64]
    Theta: NDArray[np.float64]
    Tau: NDArray[np.float64]
    Positions: NDArray[np.int64]
    Itineraries: NDArray[str]
    Labels: NDArray[np.int32]

    def get_trajectory(self, i: int) -> "Ensemble":
        """@brief Gets the `i`th row (a trajectory) of each field in `Ensemble`
        @return An `Ensemble` object with 1xM arrays
        """
        return Ensemble(**{f.name: getattr(self, f.name)[i] for f in fields(self)})

@dataclass
class PolygonalChannel:
    """@brief Holds and computes essential parameters for the billiard
    """
    d : np.float64
    alpha : np.float64
    x : NDArray[np.float64] = field(init = False)
    y : NDArray[np.float64] = field(init = False)
    edge_len : NDArray[np.float64] = field(init = False)
    cum_len : NDArray[np.float64] = field(init = False)
    def __post_init__(self):
        dx = 1/2.
        dy = dx*(1/math.tan(self.alpha/2.))
        self.x = np.array([0.,dx,2*dx,2*dx,dx,0.,0.])
        self.y = np.array([0.,dy,0.,self.d,self.d+dy,self.d,0.])
        self.edge_len = np.sqrt(np.diff(self.x)**2 + np.diff(self.y)**2)
        self.cum_len = np.insert(np.cumsum(self.edge_len), 0, 0)

def get_study_parameters(study=STUDY_PREFIX + "0"):
    """@brief 
    Used to create a `PolygonalChannel` object. See `get_polygonal_channel()` 
    """
    jconfig = json.load(open(FILE_JSON))
    return jconfig[study]["parameters"] 

def get_polygonal_channel(study=STUDY_PREFIX + "0"):
    """@brief Initialises a `PolygonalChannel` from json config
    """
    parameters = get_study_parameters(study)
    return PolygonalChannel(d=parameters["width"],alpha=parameters["alpha"])

def get_study_files(study=STUDY_PREFIX + "0", warn=False):
    """@brief Gets the list of common filenames from ../config.json
    @param warn Warn about absent files. 
    """
    jconfig = json.load(open(FILE_JSON))
    jfiles = jconfig[study]["files"] 
    fnames = jfiles.values()
    if all(Path(f).exists() for f in fnames):
        return jfiles
    else:
        missing = [f for f in fnames if not Path(f).exists()]
        if warn:
            print(f"WARNING: File(s) {missing} not found.")
        return jfiles
    #return [f for f in fnames if Path(f).exists()]

def get_ensemble(study=STUDY_PREFIX + "0"):
    """@brief Loads data into an `Ensemble` 
    @param study 
    @return Ensemble data class object
    TODO: Guard against large files?
    """
    all_files = get_study_files(study)
    H           = np.loadtxt(all_files["H"],           ndmin=2,dtype=np.float64)
    Tau         = np.loadtxt(all_files["Tau"],         ndmin=2,dtype=np.float64)
    Theta       = np.loadtxt(all_files["Theta"],       ndmin=2,dtype=np.float64)
    Positions   = np.loadtxt(all_files["Positions"],   ndmin=2,dtype=np.int64)
    Itineraries = np.loadtxt(all_files["Itineraries"], ndmin=2,dtype=str)
    Labels      = np.loadtxt(all_files["Labels"],      ndmin=2,dtype=np.int32)
    return Ensemble(H,Theta,Tau,Positions,Itineraries,Labels)

def get_trajectory_member(study=STUDY_PREFIX + "0", member_name="H"):
    """@brief Loads one of H,Tau,Theta,Positions or Itineraries 
    """
    all_files = get_study_files(study)
    ftype = np.float64
    if member_name == "Itineraries":
        ftype = str
    elif member_name == "Positions" or member_name == "Labels":
        ftype = np.int64
    return np.loadtxt(all_files[member_name],ndmin=2,dtype=ftype)

def get_ensemble_row_chunk(study, start_row, chunk_size):
    """@brief Loads a chunk (subset of what is available) into ensemble 
    @param idx Index of the starting row
    @param chunk_size Number of rows in chunk
    @param study 
    @return Ensemble data class object
    """
    all_files = get_study_files(study)
    H           = np.loadtxt(all_files["H"],           ndmin=2, skiprows=start_row, max_rows=chunk_size, dtype=np.float64)
    Tau         = np.loadtxt(all_files["Tau"],         ndmin=2, skiprows=start_row, max_rows=chunk_size, dtype=np.float64)
    Theta       = np.loadtxt(all_files["Theta"],       ndmin=2, skiprows=start_row, max_rows=chunk_size, dtype=np.float64)
    Positions   = np.loadtxt(all_files["Positions"],   ndmin=2, skiprows=start_row, max_rows=chunk_size, dtype=np.int64)
    Labels      = np.loadtxt(all_files["Labels"],      ndmin=2, skiprows=start_row, max_rows=chunk_size, dtype=np.int32)
    # BUG:
    #Itineraries = np.loadtxt(all_files["Itineraries"], ndmin=2, skiprows=start_row, max_rows=chunk_size, dtype=str)
    Itineraries = np.array(pd.read_csv(all_files["Itineraries"], skiprows=start_row, nrows=chunk_size, dtype=str, sep=" ", header=None))
    return Ensemble(H,Theta,Tau,Positions,Itineraries,Labels)


