#ifndef TRAJECTORY_H_INCLUDED
#define TRAJECTORY_H_INCLUDED

#include <vector>
#include "config.h"
#include "Particle.h"
#include "SMapHelper.h"
#include "Writer.h"
#include "json.hpp" // https://github.com/nlohmann/json
		    
/**@brief Low cost equivalent of std::vector<SParticle>
 * 
 * Implements `getItinerary` functions which are computed in post
 */
template <typename T>
struct STrajectory {

    STrajectory(const size_t& N) {
	nIterations = 0;
	H.resize(N);
	Theta.resize(N);
	Tau.resize(N);
	Position.resize(N);
	Label.resize(N);
    }

    /// Update trajectory data with an raw elements
    void Update(const T& h, const T& theta, const T& tau, const int_ll& position, const int& label) {
	H[nIterations] = h;
	Theta[nIterations] = theta;
	Tau[nIterations] = tau;
	Position[nIterations] = position;
	Label[nIterations] = label;
	nIterations++;
    }

    /// Update trajectory data with an `SParticle`
    void Update(const SParticle<T>& SP) {
	H[nIterations] = SP.H;
	Theta[nIterations] = SP.Theta;
	Tau[nIterations] = SP.Tau;
	Position[nIterations] = SP.Position;
	Label[nIterations] = SP.Label;
	nIterations++;
    }

    /**@brief 
     * @return Vector  
     */ 
    //std::vector<std::string> getItinerary() {
    //    std::vector<std::string> itin(Theta.size(),"");
    //    if (config::d > .9) {
    //        for (size_t i = 0; i < Theta.size(); ++i) 
    //    	itin[i] = getItineraryFromRegionLabel_IH(Label[i]);
    //    } else {
    //        for (size_t i = 0; i < Theta.size(); ++i) {
    //    	T theta = Theta[i]; 
    //    	if (!(Theta[i] > -PI2 && Theta[i] < PI2)) { 
    //    	    theta = PI - Theta[i]; 
    //    	}
    //    	itin[i] = getItineraryFromRegionLabel_FH<T>(H[i], theta, Label[i]);
    //        }
    //    }
    //    return itin;
    //}

    std::vector<std::string> getItinerary() {
        std::vector<std::string> itineraries(Theta.size(),"");
        if (config::widthSelection == 1) {
            for (size_t i = 0; i < Theta.size(); ++i) 
        	itineraries[i] = getItineraryFromLabel(-1., -1., Label[i]);
        } else {
            for (size_t i = 0; i < Theta.size(); ++i) {
        	T theta = Theta[i]; 
        	if (!(Theta[i] > -PI2 && Theta[i] < PI2)) { 
        	    theta = PI - Theta[i]; 
        	}
        	itineraries[i] = getItineraryFromLabel<T>(H[i], theta, Label[i]);
            }
        }
        return itineraries;
    }

    /**@brief Writes all field members
     * @param trajWriters Fixed array of writer objects 
     */ 
    void Write(std::array<Writer,6>& trajWriters) {
	trajWriters[0].WriteRowVector<T>(H);
	trajWriters[1].WriteRowVector<T>(Tau);
	trajWriters[2].WriteRowVector<T>(Theta);
	trajWriters[3].WriteRowVector<int>(Label);
	trajWriters[4].WriteRowVector<int_ll>(Position);
	trajWriters[5].WriteRowVector<std::string>(getItinerary());
    }

    std::vector<T> H, Theta, Tau;  // Entrance heights, angles and dwell times
    std::vector<int> Label;	   // Region label
    std::vector<int_ll> Position;  // Lifted position of particle
    size_t nIterations;		   // Index to each member
};

/**@brief Gets writer objects for each field member
 * Main purpose is to reduce IOPS when looping
 */ 
std::array<Writer,6> getTrajectoryWriters() {
    nlohmann::json jfiles = config::getJSONFiles();
    std::array<std::string,6> filenames = { jfiles["H"], jfiles["Tau"], jfiles["Theta"], jfiles["Labels"], 
					   jfiles["Positions"], jfiles["Itineraries"] };
    std::array<Writer,6> TrajWriters;
    for (size_t i = 0; i < TrajWriters.size(); ++i)
	TrajWriters[i] = Writer(filenames[i]);
    return TrajWriters;
}

template <typename T>
void Summarise(const STrajectory<T>& Traj) {
    //int_ll minDisplacement = - Traj.Position[0];
    //int_ll maxDisplacement = - Traj.Position[0];
    //std::vector<size_t> regionsVisited = ; 
    //T tarjectoryLength = std::cumsum(tau)
}


#endif
