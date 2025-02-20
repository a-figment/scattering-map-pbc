#ifndef WRITER_H_INCLUDED
#define WRITER_H_INCLUDED

#include <iomanip>
#include <fstream>

/*@brief Minimal wrapper class for writing data 
*/ 
class Writer
{
    private:
	std::ofstream WriteStream;
	std::string filename; // should include full path
	bool append;	      // append to file

    public:
	Writer() = default;
        Writer(const std::string& fname) : filename(fname)
        {
	    append = false;
	    WriteStream.open(filename.c_str(), std::ios::out);
	    if (!WriteStream.is_open()) 
		throw std::runtime_error("ERROR: Could not open file: " + filename + 
					 ". Check that the directory exists.");
	    //WriteStream << std::setprecision(std::numeric_limits<T>::digits10) << std::showpoint;
	    setStreamPrecision<double>();
        }

        Writer(const std::string& fname, const bool& app) : filename(fname), append(app)
        {
	    if (append)
		WriteStream.open(filename.c_str(), std::ios::out | std::ios::app);
	    else                                             
		WriteStream.open(filename.c_str(), std::ios::out);

	    if (!WriteStream.is_open()) 
		throw std::runtime_error("ERROR: Could not open file: " + filename + 
					 ". Check that the directory exists.");
	    setStreamPrecision<double>();
        }

	template <typename T>
	void setStreamPrecision() {
	    WriteStream << std::setprecision(std::numeric_limits<T>::digits10) << std::showpoint;
	}

	template <typename T>
	void WriteRowVector(const std::vector<T>& v) {
	    for (size_t i = 0; i < v.size()-1; ++i)
		WriteStream << v[i] << " ";
	    WriteStream << v.back() << std::endl;
	}

	// TODO: Handrails
	template <typename... Vectors>
	void WriteVectorsByRow(const Vectors&... v) {
	    ([&]
    	    {
	        for (const auto& elem: v) {
	            WriteStream << elem << " ";
	        }
	        WriteStream << std::endl;
    	    } (), ...);
	}

	// TODO: Handrails. Vectors must be same size
	template <typename... Vectors>
	void WriteVectorsByCol(const Vectors&... v) {
    	    const size_t size = std::get<0>(std::forward_as_tuple(v...)).size();
	    constexpr size_t pack_size = sizeof...(v);
    	    for (size_t i = 0; i < size; ++i) {
		size_t pack_idx = 0;
	        ((WriteStream << v[i] << (++pack_idx < pack_size ? " " : "")), ...);
    	        WriteStream << std::endl;
    	    }
	}
};

#endif
