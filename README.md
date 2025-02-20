
## Overview

Minimalistic C++ implementation of the scattering map introduced and defined in:

>Orchard J., Frascoli F., Rondoni., L and Mejía-Monasterio C., _Particle transport in open polygonal billiards: a scattering map_. Chaos. **2024** (34), p.123158. Available at: https://doi.org/10.1063/5.0219730.

The codebase consists of a header only implementation of the map along with helper functions/classes for writing data and generating trajectories. A collection of python scripts are also available for reconstructing billiard trajectories from scattering map data, animating 'clouds' of billiard particles and visualising the coordinate space.

![readmegif](https://github.com/user-attachments/assets/b5690e89-65c6-480e-81ee-909e62f9a3f6)

(See https://www.youtube.com/@Bills_Yard for additional animations)
## Getting started

1. Clone `git clone a-figment/scattering-map-pbc` 
2. Populate directories with `mkdir -p data/{study_0,study_1} results`
3. `cd cfiles/` and attempt `make`  
4. Run the sample code `make run TARGET=main_sample ARGS="0 0"

This will write data to `data/study_0` and generate `config.json`. Running `python ../pyfiles/sample_main.py` will generate and save some visualisations of the data (trajectory plots, animations and coordinate space) in `results/`. 


Alternatively, copy the project header files `cfiles/include/*.h` into a subdirectory of choice (e.g. `project-raw/headers`) and from `project-raw/`, populate with `mkdir -p data/{study_0,study_1} results` and compile as preferred. 

## Performance
Tabulation of raw runtime (scattering map performance only) with d = 1/2. 

| No. particle | No. Iterates | Runtime (s) | Max RSS (kb)|
|:------------:|:------------:|:-----------:|:-----------:| 
|     $10^4$   |    $10^4$    |    ~10      |     ~4032   |   
|     $10^4$   |    $10^5$    |    ~100     |     ~8112   |
|     $10^4$   |    $10^6$    |    ~999     |     ~42948  |
|              |              |             |             |
|     $10^5$   |    $10^4$    |    ~102     |     ~4423   |
|     $10^5$   |    $10^5$    |    ~1007    |     ~7778   |
|     $10^5$   |    $10^6$    |    ~10045   |     ~42884  |



