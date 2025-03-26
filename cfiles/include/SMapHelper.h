/**
 * @brief Defines and implements helper functions used throughout the project
 */

#ifndef SMAPHELPER_H_INCLUDED
#define SMAPHELPER_H_INCLUDED

#include <algorithm>
#include "config.h"
#include "CoordinateSpace.h"
#include "Random.h"

/**@brief Gets the itinerary 
 *
 * Finite horizon itineraries require (h,theta) at at runtime due to bouncing trajectories
 * @param h
 * @param theta
 * @param label Region label
 * @return string Itinerary belonging to config::itineraries
*/
template <typename T>
std::string getItineraryFromLabel(const T& h, const T& theta, const int& label) {

    if (config::widthSelection == 0) { 
    	std::string bouncing_itinerary = config::itineraries[label];
	/// Sub-regions \beta_1^(n), \beta_12^(n), \beta_18^(n)
	if (label == 1) {
    	    for (int_ll i = 0; i < alpha(config::d,h,theta); ++i)
    	        bouncing_itinerary += "04"; 
    	    return bouncing_itinerary + "R";
    	} else if (label == 12) {
    	    for (int_ll i = 0; i < kappa(config::d,h,theta); ++i)
    	        bouncing_itinerary += "31";
    	    return bouncing_itinerary + "R";
    	} else if (label == 18) {
    	    for (int_ll i = 0; i < gamma(config::d,h,theta); ++i)
    	        bouncing_itinerary += "31"; 
    	    return bouncing_itinerary + "R";
    	} else {
	    return config::itineraries[label];
	}
    } else if (config::widthSelection == 1) {
	return config::itineraries[label];
    } else {
	throw std::runtime_error("ERROR: Invalid `widthSelection`.");
    }
}

/**@brief Determines the region label corresponding to the point (h,\theta)
* @param h Entrance height
* @param theta Entrance angle in (-pi/2, pi/2)
* @return Integer label for the region. -1 signals error elsewhere
*/
template <typename T>
int whichRegion(const T& h, const T& theta) {

    T thetap = theta;
    if (!(theta > -PI2 && theta < PI2))
	thetap = PI - theta;
    // Special case for the unique region \beta_19
    if (config::widthSelection == 1) 
	if (IN_BETA19(config::d,h,thetap)) 
	    return 19; 

    if (IN_G0(config::d,h,thetap)) {
	if (config::widthSelection == 0) {
	    if (IN_BETA0(config::d,h,thetap)) { return 0; } 
	    else if (IN_BETA1(config::d,h,thetap)) { return 1; }
	    else if (IN_BETA2(config::d,h,thetap)) { return 2; }
	    else if (IN_BETA3(config::d,h,thetap)) { return 3; }
	    else if (IN_BETA4(config::d,h,thetap)) { return 4; }
	    else if (IN_BETA5(config::d,h,thetap)) { return 5; }
	    else if (IN_BETA6(config::d,h,thetap)) { return 6; }
	    else if (IN_BETA7(config::d,h,thetap)) { return 7; }
	    else if (IN_BETA8(config::d,h,thetap)) { return 8; }
	    else if (IN_BETA9(config::d,h,thetap)) { return 9; }
	    else if (IN_BETA10(config::d,h,thetap)) { return 10; }
	    else if (IN_BETA11(config::d,h,thetap)) { return 11; }
	    else { return -1; }
	} else {
	    if (IN_BETA0(config::d,h,thetap)) { return 0; } 
	    else if (IN_ZJ0(config::d,h,thetap)) { return 4; }
	    else if (IN_ZJ1(config::d,h,thetap)) { return 5; }
	    else if (IN_ZJ2(config::d,h,thetap)) { return 6; }
	    else if (IN_ZJ3(config::d,h,thetap)) { return 8; }
	    else if (IN_ZJ4(config::d,h,thetap)) { return 9; }
	    else if (IN_ZJ5(config::d,h,thetap)) { return 10; }
	    else if (IN_ZJ6(config::d,h,thetap)) { return 11; }
	    else { return -1; }
	}
    } else if (IN_G2(config::d,h,thetap)) {
	if (IN_BETA12(config::d,h,thetap)) { return 12; }
    	else if (IN_BETA13(config::d,h,thetap)) { return 13; }
    	else if (IN_BETA14(config::d,h,thetap)) { return 14; }
    	else if (IN_BETA15(config::d,h,thetap)) { return 15; }
    	else { return -1; }
    } else if (IN_G3(config::d,h,thetap)) {
	if (IN_BETA16(config::d,h,thetap)) { return 16; }
    	else if (IN_BETA17(config::d,h,thetap)) { return 17; }
    	else if (IN_BETA18(config::d,h,thetap)) { return 18; }
    	else { return -1; }
    } else { return -1; }
}


/**@brief Samples a given region via uniform rejection 
* @param label Numeric label for the region
* @param N  Number of points to get
* @param return Two-tuple of vectors, heights and angles
* TODO: Wheels spin if label is invalid
* TODO: Optimise with bounding boxes?
*/
template <typename T>
std::tuple<std::vector<T>, std::vector<T>> getPointsInRegion(const int& label, const size_t& N) {
    Random<T> h, theta;
    std::vector<T> heights(N,0);
    std::vector<T> thetas(N,0);

    // Rejection sample
    for (size_t i = 0; i < N; ++i) {
	int region_label = -1;
	T h_test = 0;
	T theta_test = 0;
	do {
	    h_test = h.getUniformRandom(0, config::d);
	    theta_test = theta.getUniformRandom(-PI/(T)2, PI/(T)2);
	    region_label = whichRegion<T>(h_test, theta_test);
	    if (region_label == -1) {
		throw std::runtime_error("ERROR: Failed to find test point in `getPointsInRegion()`.");
	    }
	} while (region_label != label);
	heights[i] = h_test;
	thetas[i] = theta_test;
    }
    return { heights, thetas };
}

#endif
