#ifndef SCATTERINGMAP_H_INCLUDED
#define SCATTERINGMAP_H_INCLUDED

#include "config.h"
#include "CoordinateSpace.h"
#include "Random.h"
#include "SMapHelper.h"
#include "Particle.h"
#include "Trajectory.h"

template <typename T>
struct ScatteringMap {
	T d;

	ScatteringMap(const T& channelWidth) : d(channelWidth) { }

	/**@brief Updates the state of a particle through an iteration of the scattering map	 
	 * 
	 * @param SP Scattering particle
	*/
	void Evolve(SParticle<T>& SP) {
	    bool dirFlag = false; // Flag for particle direction
	    SHatMap(SP, dirFlag);
	}

	/**
	 * @brief Iterates the map and records particle state at selected iterates `N`
	 *
	 * @param SP Scattering particle
	 * @param N Whole number sequence of iterates to be extracted
	 * TODO: `SP.Tau = 0`? Should be able to iterate from any starting state 
	*/
	STrajectory<T> getTrajectory(SParticle<T>& SP, const std::vector<size_t>& N) {
	    size_t j = 0;
	    STrajectory<T> STrajectory(N.size());
	    for (size_t i = 0; i < N.back(); ++i) {
		if (N[j] == i) { // Record the state of the N[j] iterate
		    STrajectory.Update(SP); 
		    j++;
		}
		Evolve(SP);
	    }
	    STrajectory.Update(SP); 
	    return STrajectory;
	}

	/**
	 * @brief Iterates the map and records particle state at selected times `Time`
	 *
	 * @param SP Scattering particle
	 * @param Time Times (lengths) at which to record the particle state
	 * @return STrajectory Scattering trajectory struct   
	 *
	 * NOTE: For d=1/2 the dwell time has no bound and the particle can be
	 * effectively localised to a particular cell for extended periods. 
	 * e.g. any given sequence of desired time points 'Time[i:j]' can be 
	 * exceeded in one call to `Evolve(...)`. In this scenario the `STrajectory` 
	 * is filled with replicates of `SP` spanning the localised time frame. 
	*/
	STrajectory<T> getTrajectory(SParticle<T>& SP, const std::vector<T> Time) {
    	
    	    // Evolve in continuous time
    	    size_t tidx = 0;  // User must specify Time[0] = 0 otherwise it gets ignored
    	    T TAU = SP.Tau;  // Current time
    	    T pTAU = TAU;    // Previous time
    	
    	    STrajectory<T> STrajectory(Time.size()); // Stores variables 
    	
    	    // Evolve map until time exceeds target time
    	    //while (TAU < Time.back()) {
    	    while (tidx < Time.size()) {
    	
		SParticle<T> pSP = SP; // Particle state prior to update
    		Evolve(SP);
    		TAU += SP.Tau;

		// Particle has jumped over at least one desired time step
		while (tidx < Time.size() && TAU >= Time[tidx]) {
		    if (abs(TAU - Time[tidx]) < abs(pTAU - Time[tidx])) {
			STrajectory.Update(SP);
		    } else { // or is previous iterate:
			STrajectory.Update(pSP);
		    }
		    tidx++;
		}
    		pTAU = TAU;
    	    }
    	    return STrajectory;
    	}

	// Eq. 1
	// @brief Iterates the scattering map 
	void SHatMap(SParticle<T>& SP, bool& dirFlag) {
	    if (SP.Theta > -PI2 && SP.Theta < PI2) { // Moving in the positive direction
		SMap(SP, dirFlag);
	    } else {
		dirFlag = true;	
		SP.Theta = PI - SP.Theta;
		SMap(SP, dirFlag);
		SP.Theta = PI - SP.Theta; // Revert 
	    }
	}

	// S(h,\theta) or S(h, \pi - \theta)
	// Note: SP.Label must always be initialised
	void SMap(SParticle<T>& SP, bool& dirFlag)  {
	    if (IN_BETA19(d,SP.H,SP.Theta)) {        // RIGHT EXIT, R
		SP.Label = 19;
		SR(SP.H, SP.Theta, SP.Tau);
	    } else if (IN_G0(d,SP.H,SP.Theta)) {     // BOTTOM LEFT, \Gamma_0
		Sg0(SP.H, SP.Theta, SP.Tau, SP.Label);
	    } else if (IN_G3(d,SP.H,SP.Theta)) {     // TOP LEFT, \Gamma_3
		Sg3(SP.H, SP.Theta, SP.Tau, SP.Label);
	    } else {					// TOP RIGHT, \Gamma_2
		Sg2(SP.H, SP.Theta, SP.Tau, SP.Label);
	    }
	    UpdatePosition(SP.Theta, SP.Position, dirFlag);
	    //UpdateLabel(SP);
	    SP.Label = whichRegion<T>(SP.H,SP.Theta);
	 }

	//void UpdateLabel(SParticle<T>& SP) {
	//    if (SP.Theta > -PI2 && SP.Theta < PI2) { // Moving in the positive direction
	//	SP.Label = whichRegion<T>(SP.H,SP.Theta);
	//    } else {
	//	SP.Label = whichRegion<T>(SP.H,PI-SP.Theta);
	//    }
	//}

	// @brief Updates the position of a particle based on its exit angle 
	void UpdatePosition(const T& theta, int_ll& x, const bool& dirFlag) {
	    if (dirFlag && (cos(theta) >= 0)) { 
		x = x - 1;
	    } else if (!dirFlag && (cos(theta) > 0)) {
		x = x + 1;
	    } else {
		x = x;
	    }
	}

	void Sg0(T& h, T& theta, T& tau, const int& currLabel) {
	    switch (currLabel) {
		case 0: {
		    tau = -h / sin(theta);
		    h = -h * (1/tan(theta));
		    theta = PI2 - theta;
		} break;
		case 1: {
		    int_ll ALPHA = alpha(d,h,theta); 
		    tau = d*(alpha(d,h,theta)+3) / cos(theta);
		    h = h + d*(ALPHA + 3)*tan(theta) + d*(ALPHA + 1);
		    theta = theta;
		} break;
		case 2: {
		    tau = -(h + d) / sin(theta);
		    h = (h+d)*(1/tan(theta)) + 5*d;
		    theta = -(theta + PI2);
		} break;
		case 3: {
		    tau = (3*DX + 2*d)/cos(theta);
		    h = h + d*(1 + 5*tan(theta));
		    theta = theta;
		} break;
		case 4: {
		    tau = (d + DX2)/cos(theta);
		    h = h + (d + DX2) * tan(theta) + d;
		    theta = theta;
		} break;
		case 5: {
		    tau = -h/sin(theta);
		    h = h*(1/tan(theta)) + 2*d + DX2;
		    theta = -(theta + PI2);
		} break;
		case 6: {
		    tau = (DX2 + 2*d)/cos(theta);
		    h = h + (DX2 + 2*d)*tan(theta);
		    theta = theta;
		} break;
		case 7: {
		    tau = (3*DX + 2*d)/cos(theta);
		    h = h + (2 + d)*tan(theta) - d;
		    theta = theta;
		} break;
		case 8: {
		    tau = (DX2 + 2*d)/cos(theta);
		    h = DX2 - (h + (DX2 + 2*d) * tan(theta));
		    theta = theta + PI;
		} break;
		case 9: {
		    tau = (DX2 + 2*d)/cos(theta);
		    h = DX2 - (h + (DX2 + 2*d)*tan(theta));
		    theta = theta + PI;
		} break;
		case 10: {
		    tau = (DX2-h)/sin(theta);
		    h = DX2 + 2*d - (DX2 - h)/tan(theta);
		    theta = 3*PI2 - theta;
		} break;
		case 11: {
		    tau = (DX2 + d)/cos(theta);
		    h = d + DX2 - (h + (DX2 + d)*tan(theta));
		    theta = theta + PI;
		} break;
		default:
		    throw std::runtime_error("Label did not correspond to any region in `Sg0(...)`.");
		    break;
	    }
	}

	// @brief Returns index of the Sg2 branch
	void Sg2(T& h, T& theta, T& tau, const int& currLabel) {
	    switch (currLabel) {
		case 12: {
		    tau = (DX2+kappa(d,h,theta)*d)/cos(theta);
		    h = (h + (DX2 + kappa(d,h,theta)*d)*tan(theta)) - d*kappa(d,h,theta);
		    theta = theta;
		} break;
		case 13: {
		    tau = (DX2 + d)/cos(theta);
		    h = d + 1 - (h + (DX2 + d)*tan(theta));
		    theta = theta + PI;
		} break;
		case 14: {
		    tau = (d + DX2 - h)/sin(theta);
		    h = (DX2 + d) - (1/tan(theta)) * (d + DX2 - h);
		    theta = -theta + 3*PI2;
		} break;
		case 15: {
		    tau = DX2/cos(theta);
		    h = 1 + 2*d - (h + tan(theta));
		    theta = theta + PI;
		} break;
		default:
		    throw std::runtime_error("Label did not correspond to any region in `Sg2(...)`.");
		    break;
	    }
	}

	// @brief Returns index of the Sg3 branch
	void Sg3(T& h, T& theta, T& tau, const int& currLabel) {
	    switch (currLabel) {
		case 16: {
		    tau = DX2/cos(theta);
		    h = 1 + 2*d - (h + tan(theta));
		    theta = theta + PI;
		} break;
		case 17: {
		    tau = (1 + 2*d - h)/sin(theta);
		    h = DX2 - (1/tan(theta))*(1 + 2*d - h);
		    theta = (3*PI2) - theta;
		} break;
		case 18: {
		    tau = (DX2 + d*(1 + gamma(d,h,theta)) - h)/sin(theta);
		    h = (1/tan(theta)) * (DX2 + d*(1 + gamma(d,h,theta)) - h) - d*(gamma(d,h,theta) - 1);
		    theta = PI2 - theta;
		} break;
		default:
		    throw std::runtime_error("Label did not correspond to any region in `Sg3(...)`.");
		    break;
	    }
	}

	// For infinite horizon, we can exit without collision
	void SR(T& h, T& theta, T& tau) {
	    tau = 1/cos(theta);
	    h = tan(theta) + h;
	    theta = theta;
	}
};

#endif
