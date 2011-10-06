#include "histogram.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <functional>

using std::cout;
using std::string;
using std::fstream;

namespace cg {
    
    // Initializes a histogram with zero bins.
    Histogram::Histogram(double x_min, double x_max, double dx) {
        double range = x_max - x_min;
        int n = (int)ceil(range/dx);
        _x0 = x_min;
        _dx = range / double(n);
        _dxi = 1.0/_dx;
        histogram.assign(n, 0.0);
    }

    //! Adds a new value to a bin in the histogram.
    void Histogram::add(double x) {
        size_t bin = (size_t)floor((x - _x0) * _dxi);
        if (bin < histogram.size()) histogram[bin] += 1.0;
    }

    // Writes a histogram to a file readable by gnuplot.
    void Histogram::write(string path) const {                
        std::fstream fid(path, std::ios::out);
        double x(_x0 + 0.5*_dx);
        for (auto h=histogram.begin(); h!=histogram.end(); ++h, x+=_dx) {
            fid << x << "\t" << *h << "\n";
        }
    }

    // Scales the histogram by a factor and by radial density.
    void Histogram::scale(double s, bool radial_density) {
        static const double four_thirds_pi = 4.0/3.0*acos(-1.0);
        double x(_x0);
        for (auto hi=histogram.begin(); hi!=histogram.end(); ++hi) {
            if (radial_density) {
                double xp = x + _dx;
                double vshell = four_thirds_pi*(xp*xp*xp - x*x*x);
                // If RDF then we want count per shell volume.
                *hi *= s/vshell;
                x = xp;
            }
            else *hi *= s;
        }
    }

    // Adds two (compatible) histograms.
    Histogram& Histogram::operator+=(const Histogram &y) {
        if (histogram.empty()) {
            *this = y;
            return *this;
        }
        if (y.histogram.size() != histogram.size()) 
            cout << "Error in histogram add: wrong # of bins.\n";
        else if (y._dx != _dx) cout << "Error in histogram add: wrong dx.\n";
        else if (y._x0 != _x0) cout << "Error in histogram add: wrong x0.\n";
        else {            
            auto hi = histogram.begin();
            auto yi = y.histogram.cbegin();
            for (; hi!=histogram.end(); ++hi, ++yi) {
                *hi += *yi;
            }
        }
        return *this;
    }
   
    //! Divides the histogram by the sum of its contents.
    void Histogram::normalize() {
        double d=std::accumulate(histogram.cbegin(),histogram.cend(),0.0);
        if (d == 0) {
            cout << "Warning: histogram is empty.\n";
            return;
        }
        else d = _dxi / d;
        for (auto h=histogram.begin(); h!=histogram.end(); ++h) *h *= d;                
    }
}
