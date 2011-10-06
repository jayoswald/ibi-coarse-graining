#pragma once
#include <vector>
#include <string>
#include <map>

namespace cg {

class LAMMPS_Data;

//! A histogram with a bin size and starting bin position.
class Histogram {
public:
    //! Initializes a histogram with zero bins.
    Histogram(double x_min=0.0, double x_max=0.0, double dx=1.0);
    //! Adds a new value to a bin in the histogram.
    void add(double x);
    //! Writes a histogram to a file readable by gnuplot.
    void write(std::string path) const;
    //! Scales the histogram by a factor and by radial density.
    void scale(double s, bool radial_density);
    //! Divides the histogram by the sum of its contents.
    void normalize();

    Histogram& operator+=(const Histogram &y);
    //! Returns the number of bins.
    size_t size() const { return histogram.size(); }

    std::vector<double> histogram;
private:
    double _x0, _dx, _dxi;    
};

}