#ifndef RANDOM_H_INCLUDED
#define RANDOM_H_INCLUDED

//#include <boost/multiprecision/cpp_dec_float.hpp>
//#include <boost/multiprecision/cpp_int.hpp>
//#include <boost/math/constants/constants.hpp>
#include <random>
#include <boost/random.hpp>

template <typename T>
class Random
{
    public:
	//using ibits = boost::random::independent_bits_engine;
	//using genny = ibits<boost::mt19937, std::numeric_limits<T>::digits, boost::multiprecision::cpp_int>;
	using genny = boost::random::mt19937;
	Random() {
    	    std::random_device rd;
    	    gen = genny(rd());
    	}
    	T getUniformRandom(const T &a, const T &b) {
    	    ur = boost::random::uniform_real_distribution<T>(a, b);
    	    return ur(gen);
    	}
    	T getNormalRandom(const T &mean, const T &var) {
    	    norm = boost::random::normal_distribution<T>(mean, var);
    	    return norm(gen);
    	}
    private:
        boost::random::uniform_real_distribution<T> ur;
        boost::random::normal_distribution<T> norm;
        genny gen;
};

#endif // RANDOM_H_INCLUDED
