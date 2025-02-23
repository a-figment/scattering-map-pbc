
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
4. Run the sample code `make run TARGET=main_sample ARGS="0 0"`

This will write data to `data/study_0` and generate `config.json`. Running `python ../pyfiles/sample_main.py` will generate and save some visualisations of the data (trajectory plots, animations and coordinate space) in `results/`. Snippet of `sample_main.cpp` demonstrates majority of currently available features:

```
    config::configure_runtime(argc, argv);     // <--- Run time configuration 
    //config::configure_compiletime(0, 0);     // <--- Compile time configuration (alpha=pi/2, d = 1/2)
    //config::configure_compiletime(0, 1);     // <--- Compile time configuration (alpha=pi/2, d = 1)

    ScatteringMap<float_> Map(config::d);      // <--- Implements Eq. 10,12,14,15, \hat{S}(h,\theta)
    //SParticle<float_> myParticle(0.25,PI/3.);// <--- Manual initialisation: (h,\theta) = (1/4,pi/3)
    //SParticle<float_> myParticle();          // <--- Random initialisation in coordinate space
    SParticle<float_> myParticle(6);           // <--- Random initialisation in region 6
    myParticle.Print();				   

    //// Evolve the map as:
    Map.Evolve(myParticle);                    // <--- Equivalent to \hat{S}(h,\theta)
    myParticle.Print();				   
    Map.Evolve(myParticle);                    // <--- \hat{S}^2(h,\theta)
    myParticle.Print();							   

    // Or, get a 'discrete' trajectory 
    myParticle = SParticle<float_>();
    STrajectory<float_> _Traj = Map.getTrajectory(myParticle, {0,5,10,15}); // <--- Record every 5 iterates {(h,\theta),S^5(h,\theta), ...}
        								     
    myParticle = SParticle<float_>(6);
    STrajectory<float_> Traj = Map.getTrajectory(myParticle, {0,1,2,3,4});  // <--- Record first 5 iterates starting from region 6 
    
    // Write selected trajectory data with wrapper class `Writer`
    Writer wHTL(config::DataPath + "HTL.dat");			   
    wHTL.WriteVectorsByRow(Traj.H, Traj.Theta, Traj.Label);		            // <--- Writes 3x5 plain text matrix of heights, angles, labels
        								    
    Writer wIP(config::DataPath + "IP.dat");
    wIP.setStreamPrecision<float>();				                        // <--- Manually specify stream precision (default: double)
    wIP.WriteVectorsByCol(Traj.getItinerary(), Traj.Position);		        // <--- Writes 5x2 data matrix. Itineraries are computed in post
        								         		     
    // Write trajectory data for particles starting in selected regions
    std::vector<int> myLabels = {};					                        // <--- Leave empty to get points in all regions, else select your own 
    std::vector<float_> myWeights = {};					                    // <--- Leave empty to get region areas as weights, else select your own 
    size_t nIterates = 5;						                            // <--- No. iterates of the map starting from 0 
    size_t nPoints = 1000;						                            // <--- No. of points per region (* region weight)
    WriteTrajectoryData<float_>(nPoints, nIterates, myLabels, myWeights);   // <--- Generate 10,000 trajectories per region 
    
    // Write trajectory data for particles starting in a rectangle
    std::pair<float_, float_> h_int(0.1, 0.4);
    std::pair<float_, float_> theta_int(-PI2 + EPSILON, PI/3.);	    
    WriteTrajectoryData<float_>(nPoints, nIterates, h_int, theta_int);      // <--- Initial conditions form a rectangle [0.1,0.4] x [-pi/2, pi/3]
    
    // Write coordinate space points 
    myLabels = {5,7,11};						    
    myWeights = {};						    
    nPoints = 10000;
    WriteRegionPoints<float_>(nPoints, myLabels, myWeights);		        // <--- Writes 10000 unweighted points belonging to each region {5,7,11}
```



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


## Cite

 @article{Orchard_2024, title={Particle transport in open polygonal billiards: A scattering map}, volume={34}, ISSN={1089-7682}, url={http://dx.doi.org/10.1063/5.0219730}, DOI={10.1063/5.0219730}, number={12}, journal={Chaos: An Interdisciplinary Journal of Nonlinear Science}, publisher={AIP Publishing}, author={Orchard, Jordan and Frascoli, Federico and Rondoni, Lamberto and Mejía-Monasterio, Carlos}, year={2024}, month=dec }
