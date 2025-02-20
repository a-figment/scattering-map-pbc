/**
 * @brief Defines and implements helper functions used for writing commonly needed data sets
 */
#ifndef PRODUCTION_H_INCLUDED
#define PRODUCTION_H_INCLUDED

#include <algorithm>
#include "config.h"
#include "Random.h"
#include "Writer.h"	   
#include "ScatteringMap.h"
#include "SMapHelper.h"
#include "Particle.h"
#include "Trajectory.h"

/**@brief Rejection samples points in a given region \beta_j, writing them to a file
 * 
 * @param nPoints Number of points to generate
 * @param myLabels User specified labels indicating which regions should be generated
*/
template <typename T>
void WriteRegionPoints(const int& nPoints, std::vector<int> myLabels, std::vector<T> myWeights) {
    if (myLabels.size() == 0) {
	myLabels = config::regionLabels;
    } else { // Labels must be a subset of config::all_labels[widthSelection]
	if (!std::ranges::includes(config::regionLabels, myLabels)) {
	    throw std::invalid_argument("ERROR: Invalid `myLabels` argument in @WriteRegionPoints().");
	}
    }
    if (myWeights.size() == 0) { myWeights = config::regionAreas; } 

    nlohmann::json jsonFiles = config::getJSONFiles();
    Writer region_heights(jsonFiles["Regions-H"]);
    Writer region_angles(jsonFiles["Regions-Theta"]);

    std::cout << "Writing " << nPoints << " points in region(s) " << std::flush; 
    for (size_t i = 0; i < myLabels.size(); ++i) {
	(i == myLabels.size()-1) ? std::cout << std::to_string(myLabels[i]) : std::cout 
					     << std::to_string(myLabels[i]) << ", " << std::flush;
	std::vector<T> heights_in_region, angles_in_region;
	size_t weightedPoints = std::round(nPoints * myWeights[myLabels[i]]);
	std::tie(heights_in_region, angles_in_region) = getPointsInRegion<T>(myLabels[i],weightedPoints);
	region_heights.WriteRowVector<T>(heights_in_region);
	region_angles.WriteRowVector<T>(angles_in_region);
    }
    std::cout << ".\nFiles written: " << jsonFiles["Regions-H"] << ", " << jsonFiles["Regions-Theta"] << "\n";
}

/**
 * @brief Writes trajectory data used for plotting and animating trajectories
 *
 * Total number of particles considered 
 * @param nPoints Number of points per region to generate
 * @param myLabels User specified labels indicating which regions should be generated
 * @param myWeights User specified weightings. Multiplies `nPoints` per region
 * TODO: 
*/
template <typename T> 
void WriteTrajectoryData(const size_t& nPoints, const size_t& nIterates, std::vector<int> myLabels, std::vector<T> myWeights) {
    if (myLabels.size() == 0) {
	myLabels = config::regionLabels;
    } else { // Labels must be a subset of config::all_labels[widthSelection]
	if (!std::ranges::includes(config::regionLabels, myLabels)) {
	    throw std::invalid_argument("Invalid myLabels argument in @WriteTrajectoryData()");
	}
    }
    if (myWeights.size() == 0) { myWeights = config::regionAreas; } 
    // Iterates start from 0
    std::vector<size_t> iterates(nIterates);
    std::iota(std::begin(iterates), std::end(iterates), 0);

    auto TrajWriters = getTrajectoryWriters(); // Opens files for all trajectory data
    size_t totalPoints = 0;
    std::cout << "Writing " << nIterates << " iterate trajectories from regions ";
    for (size_t i = 0; i < myLabels.size(); ++i) {
	(i == myLabels.size()-1) ? std::cout << std::to_string(myLabels[i]) : std::cout 
					     << std::to_string(myLabels[i]) << ", " << std::flush;
	std::vector<T> heights_in_region, angles_in_region;
	size_t weightedPoints = std::round(nPoints * myWeights[myLabels[i]]);
	totalPoints += weightedPoints;
	std::tie(heights_in_region, angles_in_region) = getPointsInRegion<T>(myLabels[i],weightedPoints);
	for (size_t j = 0; j < heights_in_region.size(); ++j) {
	    ScatteringMap<T> Map(config::d);
	    SParticle<T> Particle(heights_in_region[j],angles_in_region[j],0.,0);
	    STrajectory<T> Traj = Map.getTrajectory(Particle, iterates);     
	    Traj.Write(TrajWriters); // Writes each member of STrajectory to corresponding file
	}
    }
    std::cout << ". Done!\n" << totalPoints << " total trajectories written.\n";
}

/**
 * @brief Writes trajectory data for initial conditions in a rectangle [a,b]x[c,d] 
 *
 * Overloaded function
*/
template <typename T> 
void WriteTrajectoryData(const size_t& nPoints, const size_t& nIterates, 
			 const std::pair<T,T>& h_interval, const std::pair<T,T>& theta_interval) {

    if (h_interval.first < (T)0 || h_interval.second > config::d) {
	throw std::invalid_argument("Invalid height interval in @WriteTrajectoryData()");
    }
    // Iterates start from 0
    std::vector<size_t> iterates(nIterates);
    std::iota(std::begin(iterates), std::end(iterates), 0);

    auto TrajWriters = getTrajectoryWriters(); // Opens files for all trajectory data

    std::cout << "Writing " << std::to_string(nPoints) << " trajectories (" << nIterates 
	      << " iterates) from the region " << "[" << h_interval.first << ", " 
	      << h_interval.second << "] x [" << theta_interval.first << ", " << theta_interval.second << "]";
    // Generate initial condition in rectangle, produce trajectory and write
    for (size_t i = 0; i < nPoints; ++i) {
	ScatteringMap<T> Map(config::d);
	SParticle<T> P = SParticle(h_interval, theta_interval);
	STrajectory<T> Traj = Map.getTrajectory(P, iterates);     
	Traj.Write(TrajWriters); // Writes each member of STrajectory to corresponding file
    }
    std::cout << ". Done!\n";
}


#endif
