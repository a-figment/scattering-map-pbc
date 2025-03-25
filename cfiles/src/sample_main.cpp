#include <iostream>
#include <chrono>
#include <cmath>
#include <vector>

#include "json.hpp"	   // https://github.com/nlohmann/json
		    
#include "config.h"        // Handles param. selection and files
#include "Random.h"        // Random num. gen. wrapper class
#include "Writer.h"	   // ofstream wrapper class 
#include "SMapHelper.h"    // Misc. helper functions
#include "Production.h"    // Misc. funcs for writing common data sets (optional)
#include "ScatteringMap.h" // The map

int main(int argc, char *argv[])
{
    using clock = std::chrono::high_resolution_clock;
    clock::time_point t0 = clock::now();
    config::configure_runtime(argc, argv);	   // <--- Run time configuration 
    //config::configure_compiletime(0, 0);	   // <--- Compile time configuration (alpha=pi/2, d = 1/2)
    //config::configure_compiletime(0, 1);	   // <--- Compile time configuration (alpha=pi/2, d = 1)

    // Example usage below - 
    ScatteringMap<float_> Map(config::d);	   // <--- Implements Eq. 10,12,14,15, \hat{S}(h,\theta)
    //SParticle<float_> myParticle(0.25,PI/3.);	   // <--- Manual initialisation: (h,\theta) = (1/4,pi/3)
    //SParticle<float_> myParticle();		   // <--- Random initialisation in coordinate space
    SParticle<float_> myParticle(6);		   // <--- Random initialisation in region 6
    myParticle.Print();				   // <--- Write current state to console

    //// Evolve the map as:
    Map.Evolve(myParticle);			   // <--- Equivalent to \hat{S}(h,\theta)
    myParticle.Print();				   
    Map.Evolve(myParticle); 			   // <--- \hat{S}^2(h,\theta)
    myParticle.Print();							   

    // Or, get a 'discrete' trajectory 
    myParticle = SParticle<float_>();
    std::vector<size_t> iterates = {0,5,10,15};
    STrajectory<float_> _Traj = Map.getTrajectory(myParticle, iterates); // <--- Record every 5 iterates {(h,\theta),S^5(h,\theta), ...}
        								     
    myParticle = SParticle<float_>(6);
    std::vector<float_> times = {0.,10.,35.};
    STrajectory<float_> Traj = Map.getTrajectory(myParticle, times);  // <--- Record first 5 iterates starting from region 6 
    
    // Write selected trajectory data with wrapper class `Writer`
    Writer wHTL(config::DataPath + "HTL.dat");			   
    wHTL.WriteVectorsByRow(Traj.H, Traj.Theta, Traj.Label);		    // <--- Writes 3x5 plain text matrix of heights, angles, labels
        								    
    Writer wIP(config::DataPath + "IP.dat");
    wIP.setStreamPrecision<float>();				            // <--- Manually specify stream precision (default: double)
    wIP.WriteVectorsByCol(Traj.getItinerary(), Traj.Position);		    // <--- Writes 5x2 data matrix. Itineraries are computed in post
        								         		     
    // Write trajectory data for particles starting in selected regions
    std::vector<int> myLabels = {};					    // <--- Leave empty to get points in all regions, else select your own 
    std::vector<float_> myWeights = {};					    // <--- Leave empty to get region areas as weights, else select your own 
    size_t nIterates = 5;						    // <--- No. iterates of the map starting from 0 
    size_t nPoints = 1000;						    // <--- No. of points per region (* region weight)
    WriteTrajectoryData<float_>(nPoints, nIterates, myLabels, myWeights);   // <--- Generate 10,000 trajectories per region 
    
    // Write trajectory data for particles starting in a rectangle
    std::pair<float_, float_> h_int(0.1, 0.4);
    std::pair<float_, float_> theta_int(-PI2 + EPSILON, PI/3.);	    
    WriteTrajectoryData<float_>(nPoints, nIterates, h_int, theta_int);      // <--- Initial conditions form a rectangle [0.1,0.4] x [-pi/2, pi/3]
    
    // Write coordinate space points 
    myLabels = {5,9,11};						    
    myWeights = {};						    
    nPoints = 10000;
    WriteRegionPoints<float_>(nPoints, myLabels, myWeights);		    // <--- Writes 10000 unweighted points belonging to each region {5,7,11}


    clock::time_point t1 = clock::now();
    const std::string RunTime = "\033[1;32mRun duration: " 
			      + std::to_string(std::chrono::duration_cast<std::chrono::minutes>(t1 - t0).count()) 
			      + " minutes = " 
			      + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()) 
			      + " milliseconds. \n\033[0m" + "\n";
    std::cout << RunTime;
}
