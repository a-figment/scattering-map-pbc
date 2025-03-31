#ifndef TRAJECTORY_H_INCLUDED
#define TRAJECTORY_H_INCLUDED

#include <vector>
#include "config.h"
#include "Particle.h"
#include "SMapHelper.h"
#include "Writer.h"
#include "json.hpp" // https://github.com/nlohmann/json

std::string format_itinerary(const std::string& s, size_t width) {
    if (s.length() <= width) 
	return s;
    size_t keep = (width - 3) / 2;
    return s.substr(0, keep) + "..." + s.substr(s.length() - keep);
}

// TODO: Move somewhere appropriate or use boost/format, fmt
template <typename T>
void PrintVector(const std::string& name, const std::vector<T>& v, const size_t& nToPrint) {
    const size_t width = 10;
    std::cout << "    " << std::left << std::setw(width) << name << ": ";
    for (size_t i = 0; i < nToPrint && i < v.size(); ++i) 
        std::cout << std::setw(width) << v[i] << " ";
    if ((v.size() - nToPrint) > nToPrint) { std::cout << "... "; }
    for (size_t i = (v.size() > nToPrint ? v.size() - nToPrint : 0); i < v.size(); ++i) 
        if (i >= nToPrint)
	    std::cout << std::setw(width) << v[i] << " ";
    std::cout << "\n";
}  
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
	if (nIterations < Time.size() && nIterations > 0) //?
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
