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
 * TODO: Currently overspecified. 
 */
template <typename T>
struct STrajectory {

    STrajectory(const size_t& n) {
	N = n;
	nIterations = 0;
	H.resize(N);
	Theta.resize(N);
	Tau.resize(N);
	Position.resize(N);
	Label.resize(N,-1);
	Time.resize(N,(T)0);
    }    

    /// Update trajectory data with an raw elements
    void Update(const T& h, const T& theta, const T& tau, const int_ll& position, const int& label) {
	H[nIterations] = h;
	Theta[nIterations] = theta;
	Tau[nIterations] = tau;
	Position[nIterations] = position;
	Label[nIterations] = label;
	UpdateTime(tau);
	nIterations++;
    }

    /// Update trajectory data with an `SParticle`
    void Update(const SParticle<T>& SP) {
	H[nIterations] = SP.H;
	Theta[nIterations] = SP.Theta;
	Tau[nIterations] = SP.Tau;
	Position[nIterations] = SP.Position;  
	Label[nIterations] = SP.Label;
	UpdateTime(SP.Tau);
	nIterations++;
    }

    void Update(const SParticle<T>& SP, const T& cumulTime) {
	H[nIterations] = SP.H;
	Theta[nIterations] = SP.Theta;
	Tau[nIterations] = SP.Tau;
	Position[nIterations] = SP.Position;  
	Label[nIterations] = SP.Label;
	Time[nIterations] = cumulTime;
	nIterations++;
    }

    // Temporary
    void UpdateTime(const T& tau) {
	if (nIterations > 0 && nIterations < Time.size()) //?
	    Time[nIterations] = Time[nIterations-1] + tau;
    }

    /**@brief Gets currently available itineraries  
     */
    std::vector<std::string> getItinerary() const {
        std::vector<std::string> itineraries(N,"NULL");
        for (size_t i = 0; i < nIterations; ++i) 
	    itineraries[i] = getItineraryFromLabel<T>(H[i], Theta[i], Label[i]);
        return itineraries;
    }

    /**@brief Determines the lifted positions from a coordinate space itineray
     * @param labels Ordered sequence of labels / coordinate space itineray
    */
    std::vector<int_ll> getLiftedPositionsFromLabel(const std::vector<int>& labels) const {
	std::vector<int_ll> x(labels.size(),0);     // Change size to N later, fill at selected
	//std::vector<int> vx(labels.size(),1);     // Change size to N later, fill at selected
	int vx0 = 1;	    // = sgn(cos(Theta[0])) // assume initially move in positive direction
	x[0] = 0;
	int vx = vx0;
	for (size_t i = 1; i < labels.size(); ++i) {
	    std::string LorR = config::lastChars[labels[i-1]];  // Current region is trans. or refl.
	    if (LorR == "R") {
		x[i] = x[i-1] + vx;
	    } else if (LorR == "L") {
		x[i] = x[i-1];
		vx = -vx;
	    } else {
		std::cout << "ERROR: getLiftedPositionsFromLabel" << std::endl;
	    }
	}
	return x;
    }

    /**@brief Writes all field members
     * @param trajWriters Fixed array of writer objects 
     */ 
    void Write(std::array<Writer,6>& trajWriters) const {
	trajWriters[0].WriteRowVector<T>(H);
	trajWriters[1].WriteRowVector<T>(Time);
	trajWriters[2].WriteRowVector<T>(Theta);
	trajWriters[3].WriteRowVector<int>(Label);
	trajWriters[4].WriteRowVector<int_ll>(Position);
	trajWriters[5].WriteRowVector<std::string>(getItinerary());
    }    

    /**@brief Writes selected field members
     * @param trajWriters Fixed array of writer objects 
     */ 
    void WriteSelected(std::array<Writer,6>& trajWriters, const std::vector<size_t>& idx) const {
	//idx.erase(std::unique(idx.begin(), idx.end()),idx.end());
	//
	for (size_t i = 0; i < idx.size(); ++i) {
	    if (idx[i] == 0) { trajWriters[0].WriteRowVector<T>(H); }
	    else if (idx[i] == 1) { trajWriters[1].WriteRowVector<T>(Time); }
	    else if (idx[i] == 2) { trajWriters[2].WriteRowVector<T>(Theta); }
	    else if (idx[i] == 3) { trajWriters[3].WriteRowVector<int>(Label); }
	    else if (idx[i] == 4) { trajWriters[4].WriteRowVector<int_ll>(Position); }
	    else { trajWriters[5].WriteRowVector<std::string>(getItinerary()); }
	}
    }

    /**@brief Summarises 
     * @param Traj An STrajectory
     * @param nToPrint Number of entries to print in trajectory vectors
     */ 
    void Print(const size_t& nToPrint) const {
        std::vector<size_t> iterates(H.size());
        std::iota(std::begin(iterates), std::end(iterates), 0);
        std::cout << std::setprecision(2) << std::fixed << "STrajectory (d = " << config::d << ") {\n";
        std::cout << std::setprecision(4) << std::fixed;
        PrintVector("i", iterates, nToPrint);
        PrintVector<T>("H", H, nToPrint);
        PrintVector<T>("Theta", Theta, nToPrint);
        PrintVector<T>("Tau", Tau, nToPrint);
        PrintVector<T>("Time", Time, nToPrint);
        PrintVector<int_ll>("Positions", Position, nToPrint);
        PrintVector<int>("Label", Label, nToPrint);
        PrintVector<std::string>("Itinerary", getItinerary(), nToPrint);
        std::cout << "}\n";

	// Trajectory span
    	auto mM = std::minmax_element(Position.begin(), Position.end());
    	int midx = std::distance(Position.begin(), mM.first);
    	int Midx = std::distance(Position.begin(), mM.second);
	std::cout << "    Smallest displacement of x = " << *mM.first  << " occurred at time t ~ " << Time[midx] << std::endl;
	std::cout << "    Greatest displacement of x = " << *mM.second << " occurred at time t ~ " << Time[Midx] << std::endl;
    
	std::map<int, int> bin;
    	for (int label : Label) 
    	    bin[label]++;
	std::cout << "    Number of visitations to each label [Label : Count]:" << std::endl;
	std::cout << "        ";
	for (const auto& [label, count] : bin) 
    	    std::cout << "[ " << label << " : " << count << " ] ";
	std::cout << "\n";
    }

    std::vector<T> H, Theta, Tau;  // Entrance heights, angles and dwell times
    std::vector<T> Time;	   // Running sum of dwell times
    std::vector<int> Label;	   // Region label
    std::vector<int_ll> Position;  // Lifted position of particle
    size_t nIterations;		   // Index to each member (current size)
    size_t N;			   // Final size
};

/**@brief Gets writer objects for each field member
 * Main purpose is to reduce IOPS when looping
 */ 
std::array<Writer,6> getTrajectoryWriters() {
    nlohmann::json jfiles = config::getJSONFiles();
    std::array<std::string,6> filenames = { jfiles["H"], jfiles["Time"], jfiles["Theta"], jfiles["Labels"], 
					   jfiles["Positions"], jfiles["Itineraries"] };
    std::array<Writer,6> TrajWriters;
    for (size_t i = 0; i < TrajWriters.size(); ++i)
	TrajWriters[i] = Writer(filenames[i]);
    return TrajWriters;
}






#endif
