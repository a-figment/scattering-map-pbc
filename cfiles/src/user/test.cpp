
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
    return 0;
}
