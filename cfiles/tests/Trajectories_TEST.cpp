/*@brief Generic tests of trajectories
*/ 

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
#define BOOST_TEST_MODULE Trajectories
#include <boost/test/unit_test.hpp>


/**@brief Checks for stability of some known periodic trajectories
 * PERIODIC TRAJ 1: (h, {0,pi}) with 0<h<\delta y
 */
BOOST_AUTO_TEST_CASE(ballistic_traj_1) {
    std::vector<size_t> iterates = { 0, 100, 10000, 1000000, 10000000 };
    for (size_t i = 0; i < 2; ++i) {
	config::configure_compiletime(0, i);
    	ScatteringMap<float_> SM(config::d);

    	// Initial conditions
    	float_ h_random = runif<float_>(0., 1/(float_)2); 
    	std::vector<float_> theta_periodic = { 0. , PI };

    	// Particle and trajectory generation
	for (size_t j = 0; j < theta_periodic.size(); ++j) {
	   SParticle<float_> SP1(h_random,theta_periodic[j]);
    	   STrajectory<float_> ballistic_PO = SM.getTrajectory(SP1, iterates); // Horizontal periodic traj
    	   std::vector<std::string> itin = ballistic_PO.getItinerary();
	   std::cout << "Testing P.O. with initial condition h = " << h_random << ", \\theta = " << theta_periodic[j] << std::endl;

    	   for (size_t k = 0; k < iterates.size(); ++k) {
    	       BOOST_TEST((ballistic_PO.Position[k] == (int_ll)iterates[k] || ballistic_PO.Position[k] == -(int_ll)iterates[k]));
    	       BOOST_TEST((ballistic_PO.Theta[k] == 0. || ballistic_PO.Theta[k] == PI), boost::test_tools::tolerance(1e-6));
    	       BOOST_TEST(itin[k] == "L0431R");
    	       BOOST_TEST((ballistic_PO.Time[k] == (1+2*config::d)*((float_)iterates[k])), boost::test_tools::tolerance(1e-6));
    	   }
	}
    }
}

// PERIODIC TRAJ 2: (h, {pi/4,3*pi/4}) with 0<h<\delta y
BOOST_AUTO_TEST_CASE(localised_traj_1) {
    std::vector<size_t> iterates = { 0, 100, 10000, 1000000, 10000000 };
    for (size_t i = 0; i < 2; ++i) {
        config::configure_compiletime(0, i);
        ScatteringMap<float_> SM(config::d);
    
        float_ h_random = runif(0., config::d); 
        std::vector<float_> theta_periodic = { PI/(float_)4 , 3*PI/(float_)4 };

	for (size_t j = 0; j < theta_periodic.size(); ++j) {
	    std::cout << "Testing P.O. with initial condition h = " << h_random << ", \\theta = " << theta_periodic[j] << std::endl;
            SParticle<float_> SP2(h_random,theta_periodic[j]);
            STrajectory<float_> localised_PO = SM.getTrajectory(SP2, iterates); 
            std::vector<std::string> itin = localised_PO.getItinerary();
            for (size_t k = 0; k < iterates.size(); ++k) {
		BOOST_TEST(localised_PO.Position[k] == 0);
    	    	BOOST_TEST(localised_PO.H[k] == h_random, boost::test_tools::tolerance(1e-6));
		BOOST_TEST(itin[i] == "L3L");
            }
	}
    }
}

// PERIODIC TRAJ 3: (h, {-pi/4,3*pi/4}) with 0<h<\delta y
BOOST_AUTO_TEST_CASE(localised_traj_2) {
    std::vector<size_t> iterates = { 0, 100, 10000, 1000000, 10000000 };
    for (size_t i = 0; i < 2; ++i) {
        config::configure_compiletime(0, i);
        ScatteringMap<float_> SM(config::d);
    
        float_ h_random = runif(0., config::d); 
        std::vector<float_> theta_periodic = { -PI/(float_)4 , 5*PI/(float_)4};

	for (size_t j = 0; j < theta_periodic.size(); ++j) {
	    std::cout << "Testing P.O. with initial condition h = " << h_random << ", \\theta = " << theta_periodic[j] << std::endl;
            SParticle<float_> PO_particle(h_random,theta_periodic[j]);
            STrajectory<float_> localised_PO = SM.getTrajectory(PO_particle, iterates); 
            std::vector<std::string> itin = localised_PO.getItinerary();
            for (size_t k = 0; k < iterates.size(); ++k) {
		BOOST_TEST(localised_PO.Position[k] == 0);
    	    	BOOST_TEST(localised_PO.H[k] == h_random, boost::test_tools::tolerance(1e-6));
		BOOST_TEST(itin[i] == "L0L");
            }
	}
    }
}

