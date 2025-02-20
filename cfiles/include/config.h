/**
 * @brief Defines geometric constants and types, and runtime vars in a namespace 
 *
 * Required throughout the whole project. 
 */

#ifndef CONFIG_H_INCLUDED
#define CONFIG_H_INCLUDED

#include <boost/math/constants/constants.hpp>
#include <vector>
#include <sstream>
#include <fstream>
#include <filesystem>
#include "json.hpp" // https://github.com/nlohmann/json
//#include <boost/multiprecision/cpp_dec_float.hpp>
//#include <boost/multiprecision/cpp_int.hpp>
//typedef boost::multiprecision::number<mp::cpp_dec_float<50>, mp::et_off> float_50; 

//typedef float_50 float_;
//using float_ = double;
//using int_ll = std::int32_t;
//using cid_int = std::uint32_t;
typedef double float_;
typedef std::int32_t int_ll;
typedef std::uint32_t cid_int;

/// Ubiquitous consts
constexpr float_ EPSILON = std::numeric_limits<float_>::epsilon();
constexpr float_ PI = boost::math::constants::pi<float_>();
constexpr float_ PI2 = boost::math::constants::half_pi<float_>();
/// Geometric constants - used in ScatteringMap.h and CoordinateSpace.h
constexpr float_ DELTA = boost::math::constants::half<float_>(); 
constexpr float_ DX = DELTA; // \delta_x 
constexpr float_ DX2 = 1; 

namespace config {

    // Geometric and cspace parameters
    constexpr std::array<float_, 2> channel_widths = { 1/(float_)2, (float_)1 };
    constexpr std::array<float_, 1> channel_angles = { PI2 }; 
    
    // Itineraries as in tables X, Y, Z. TODO: refactor coordinatespace.h
    const std::vector<std::vector<std::string>> all_itineraries = {{ "L0L", "L04", "L04043R", "L040431R", "L04R", "L043R", 
					"L0431R", "L043131R", "L0430L", "L0340L", "L034L", "L03L", 
					"L", "L30L", "L3L", "L34L", "L43L", "L430L", "L4" },
					{ "L0L", "", "", "", "L04R", "L043R", "L0431R", "", "L0430L", 
					  "L0340L", "L034L", "L03L", "L31R", "L30L", "L3L", "L34L", 
					  "L43L", "L430L", "L431L", "LR" }};
    // Labelling convention for regions. TODO: refactor coordinatespace.h
    const std::vector<std::vector<int>> all_labels = {{ 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18 },
						      { 0,4,5,6,8,9,10,11,12,13,14,15,16,17,18,19 }};

    // Approx. area for each region. TODO: refactor coordinatespace.h
    const std::vector<std::vector<float_>> all_areas = {{0.565986, 0.0607631,  0.018407,  0.0199902, 0.0383971, 0.0199902, 
							 0.123729, 0.0199902,  0.0383971, 0.0383971, 0.0333078, 0.0274553, 
							 0.0991603,0.0274553, 0.100011, 0.0393381,  0.0393381, 0.0684178, 
							 0.192265 }, {1.13197,  -1., -1., -1.,  0.137557, 0.032183, 0.0993438, 
							 -1., 0.032183, 0.032183, 0.0446113, 0.0607631, 0.137557, 0.0607631, 
							 0.333609, 0.107756, 0.107756, 0.134789, 0.249741, 0.438825 }};
    // Directory names for output
    const std::string BaseDataPath = "../data/";
    const std::string STUDY_PREFIX = "study_"; // naming convention for data organisation
    // Common data file names 
    const std::string FILE_JSON	       = "../config.json";
    const std::string FILE_H	        = "H.dat";
    const std::string FILE_TAU	        = "Tau.dat";
    const std::string FILE_THETA        = "Theta.dat";
    const std::string FILE_LABELS       = "Labels.dat";
    const std::string FILE_POSITIONS    = "Positions.dat";
    const std::string FILE_ITINERARIES  = "Itineraries.dat";
    const std::string FILE_REGION_H     = "Regions-H.dat";     // Hoizontal cspace coords
    const std::string FILE_REGION_THETA = "Regions-Theta.dat"; // Vertical cspace coords
    
    /// ================= CONFIGURABLE VARIABLES ========================
    size_t widthSelection = 0; /// Default parameter selection is d = 1/2 
    size_t angleSelection = 0; /// and alpha = pi/2
    size_t studyID = widthSelection*channel_angles.size()+angleSelection;
    std::vector<std::string> itineraries = all_itineraries[widthSelection];
    std::vector<int> regionLabels = all_labels[widthSelection];
    std::vector<float_> regionAreas = all_areas[widthSelection];
    float_ d = channel_widths[widthSelection];
    float_ alpha = channel_angles[angleSelection];
    std::string Study = STUDY_PREFIX + std::to_string(studyID); 
    std::string DataPath = BaseDataPath + Study + "/";	
    /// ==================================================================
    
    void writeJSONConfig();

    /**@brief Initialises all configurable variables within the config::
     * @param widthSelection Indexes the available widths
     * @param angleSelection "                   " angles
    */ 
    void initialise(const size_t& angleSelect, const size_t& widthSelect) {
	if (angleSelect > channel_angles.size()-1) 
	    throw std::invalid_argument("ERROR: Aperture angle selection must be within the range [0," 
					+ std::to_string(channel_angles.size()-1) + "].");
	if (widthSelect > channel_widths.size()-1) 
	    throw std::invalid_argument("ERROR: Channel width selection must be within the range [0," 
					+ std::to_string(channel_widths.size()-1) + "].");
	widthSelection = widthSelect;
	angleSelection = angleSelect;
	d = channel_widths[widthSelection];
	alpha = channel_angles[angleSelection];
	regionLabels = all_labels[widthSelection];
	itineraries = all_itineraries[widthSelection];
	regionAreas = all_areas[widthSelection];
	
	// Study id is an integer parameterising all (alpha,d) combinations
	studyID = widthSelection*channel_angles.size()+angleSelection;
	Study = STUDY_PREFIX + std::to_string(studyID);
    	DataPath = BaseDataPath + Study + "/";
	
	writeJSONConfig();
    }

    /**@brief Generates and writes the JSON configuration file `config.json`
     * 
     * `config.json`'s purpose is to link parameter selection with a corresponding
     * list of data files.
     * 
     * NOTE: Though always written on initialise(), `config.json` contains references to 
     * all possible parameter combinations (d,alpha), regardless of runtime parameter selection
    */ 
    void writeJSONConfig() {
	std::array<std::string,8> fileNames = { FILE_H, FILE_TAU, FILE_THETA, FILE_LABELS, 
					       FILE_POSITIONS, FILE_ITINERARIES, FILE_REGION_H,
					       FILE_REGION_THETA };

	nlohmann::ordered_json jsonSTUDY;
	for (size_t i = 0; i < channel_widths.size(); ++i) {
	    for (size_t j = 0; j < channel_angles.size(); ++j) {
		std::string study = STUDY_PREFIX + std::to_string(i*channel_angles.size()+j);
	    	jsonSTUDY[study]["parameters"]["alpha"] = channel_angles[j];
	    	jsonSTUDY[study]["parameters"]["width"] = channel_widths[i];
		//jsonstudy[study]["parameters"]["region_labels"] = nlohmann::json(all_labels[i]).dump();
		jsonSTUDY[study]["parameters"]["region_labels"] = all_labels[i];
		for (size_t k = 0; k < fileNames.size(); ++k) {
		    std::string str = fileNames[k];
		    str = str.erase(fileNames[k].length() - 4);
		    jsonSTUDY[study]["files"][str] = BaseDataPath + study + "/" + fileNames[k];
		}
	    }
	}
	std::cout << "Writing ../config.json.\n";
	std::ofstream JSONWriter(FILE_JSON);
	JSONWriter << jsonSTUDY.dump(4); 
    }

    /**@brief Configuration for runtime args
     * @param widthSelection Index to the available widths
     * @param angleSelection "                   "  angles
    */ 
    void configure_runtime(int argc, char *argv[]) {
	if (argc != 3) 
	    throw std::invalid_argument("ERROR: Expected two arguments.");
	try {
	    std::stoi(argv[1]);
	    std::stoi(argv[2]);
	} catch (const std::exception& e) {
	    std::cerr << "ERROR: Invalid arguments. Expected two integers but received: " 
		<< argv[1] << " " << argv[2] << ".\n";
	}
	initialise(std::stoi(argv[1]), std::stoi(argv[2]));
    }    
    void configure_compiletime(const size_t& angleSelect, const size_t& widthSelect) {
	initialise(angleSelect, widthSelect);
    }

    /**@brief Gets files list from json object for the currently running study 
     * @return A nlohmann::json object
    */ 
    nlohmann::json getJSONFiles() {
	std::ifstream ifs(FILE_JSON);
    	nlohmann::json j = nlohmann::json::parse(ifs);
	return j[Study]["files"];
    }
}

#endif
