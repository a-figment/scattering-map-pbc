#include <iostream>
#include <chrono>
#include <cmath>
#include <vector>

#include "json.hpp"	   // https://github.com/nlohmann/json
#include "config.h"        
#include "Random.h"        
#include "Writer.h"	   
#include "SMapHelper.h"    
#include "Production.h"    
#include "ScatteringMap.h" 

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE SMapHelper
#include <boost/test/unit_test.hpp>

/**@brief Region labels must map correctly to region itineraries
 * TODO: Expected failure for bouncing trajectories. 
 */
BOOST_AUTO_TEST_CASE(itinerary_test) {
    
    // Loop over each configuration
    for (size_t i = 0; i < config::channel_widths.size(); ++i) {
	for (size_t j = 0; j < config::channel_angles.size(); ++j) {
	    config::configure_compiletime(j, i);
	    for (auto label : config::regionLabels) {
		auto [v1, v2] = getPointsInRegion<double>(label, 1);
		double h = v1[0];
		double theta = v2[0];
		std::string test_itinerary = getItineraryFromLabel(h, theta, label);
	    	BOOST_TEST(test_itinerary == config::itineraries[label]);
	    }
	}
    }
}
