#ifndef PARTICLE_H_INCLUDED
#define PARTICLE_H_INCLUDED

#include <iostream>
#include <vector>
#include "config.h"
#include "Random.h"
#include "SMapHelper.h"

/**@brief Holds information pertaining to the state of a particle
 * This is the object acted on by the scattering map
 *
 * TODO: Include dwell time to initial state after CSpace refactor 
 */

template <typename T> 
struct SParticle {
    // Uniform random state
    SParticle() {
	Random<T> h, theta;
	T H = h.getUniformRandom(0, d);
	T Theta = theta.getUniformRandom(-PI/(T)2, PI/(T)2);
	Set(H, Theta, (T)0, 0, whichRegion<T>(H,Theta));
    }

    // Hand selected particle state
    SParticle(const T& h, const T& theta, const T& tau, const int_ll& position) {
	Set(h, theta, tau, position, whichRegion<T>(h,theta));
    }

    SParticle(const T& h, const T& theta) {
	Set(h, theta, (T)0, 0, whichRegion<T>(h,theta));
    }

    SParticle(const std::pair<T,T>& h_interval, const std::pair<T,T>& theta_interval) {
	Random<T> h, theta;
	T H = h.getUniformRandom(h_interval.first, h_interval.second);
	T Theta = theta.getUniformRandom(theta_interval.first, theta_interval.second);
	Set(H, Theta, (T)0, 0, whichRegion<T>(H,Theta));
    }

    // Random inside a region 
    SParticle(const int& label) {
	auto iter = find(config::regionLabels.begin(), config::regionLabels.end(), label);
	if (iter != config::regionLabels.end()) { // Valid label
	    std::vector<T> h, theta;
	    std::tie(h, theta) = getPointsInRegion<T>(label,1);
	    Set(h[0], theta[0], (T)0, 0, label);
	} else {
	    throw std::invalid_argument("ERROR: Invalid `label` argument in @SParticle().");
	}
    }

    void Set(const T& h, const T& theta, const T& tau, const int_ll& position, const int& label) {
	H = h;
	Theta = theta;
	Tau = tau;
	Position = position;
	Label = label;
    }

    void Print() {
	std::cout << std::setprecision(3) << std::fixed
	<< "SParticle (d = " << d << ") {"
	<< "\n    H         : " << H
	<< "\n    Theta     : " << Theta << " (" << Theta*180/3.14 << " deg)"
	<< "\n    Tau       : " << Tau 
	<< "\n    Position  : " << Position
	<< "\n    Label     : " << Label
	<< "\n    Itinerary : " << getItineraryFromLabel<T>(H, Theta, Label)
	<< "\n}\n";
    }

    T d = config::d;
    T H, Theta, Tau;
    int_ll Position;
    int Label;
};

#endif
