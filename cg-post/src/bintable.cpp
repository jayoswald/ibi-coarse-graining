#include "bintable.h"
#include "lammps_data.h"
#include "cgsystem.h"
#include <cmath>
#include <set>
#include <iostream>

namespace cg {

    BinTable::BinTable(double cutoff, int step, const CgSystem *cg) { 
        // Keep pointer reference to cgsystem so we can look up bins later.
        _cgsys       = cg;
        _cutoff      = cutoff;
        _step        = step;
        _granularity = 5;

        // Number of bins along each direction.
        // Granulity to make bins smaller.
        _xbins = int(cg->box_dx(step)/cutoff);
        _ybins = int(cg->box_dy(step)/cutoff);
        _zbins = int(cg->box_dz(step)/cutoff);
        _xbins = std::max(_xbins, 1) * _granularity;
        _ybins = std::max(_ybins, 1) * _granularity;
        _zbins = std::max(_zbins, 1) * _granularity;
        
        _bins.assign(_xbins*_ybins*_zbins, std::vector<int>());

        for (int i=0; i< cg->num_beads(); ++i) {
            _bins[_find_bin(i)].push_back(i);
        }
    }

    // Returns an array of all atoms that can be within cutoff of i.
    // This might not be the best implimentation because a lot of 
    // vectors need to be copied.
    std::vector<const std::vector<int>*> BinTable::neighbors(int i) const {
        // Scaled coordinate of the atom.

        int bi  = _find_bin(i);
        std::vector<const std::vector<int>*> bins;

        int bx = bi % _xbins;
        int by = (bi/_xbins)%_ybins;
        int bz = (bi/_xbins/_ybins);

        int N = _granularity;
        for (int i=-N; i<=N; ++i) {
            int nx = (bx+i+_xbins) % _xbins;
            // Prevent wrap around.
            if (i>0 &&  nx == bx-N+_xbins) continue; 
        
            for (int j=-N; j<=N; ++j) {                
                int ny = (by+j+_ybins) % _ybins;
                if (j>0 && ny == by-N+_ybins) continue;
                for (int k=-N; k<=N; ++k) {
                    
                    int nz = (bz+k+_zbins) % _zbins;
                    if (k>0 &&nz == bz-N+_zbins) continue;
                    int nb = nx + ny*_xbins + nz*_xbins*_ybins;
                    bins.push_back(&_bins[nb]);
                }
            }
        }
        return bins; 
    }

    // Returns the index of the bin.
    int BinTable::_find_bin(int i) const {
        double dx = _cgsys->box_dx(_step);
        double dy = _cgsys->box_dy(_step);
        double dz = _cgsys->box_dz(_step);

        auto r = _cgsys->bead_position(_step, i);
        int bx  = int(double(_xbins)*r.x/dx) % _xbins;
        int by  = int(double(_ybins)*r.y/dy) % _ybins;
        int bz  = int(double(_zbins)*r.z/dz) % _zbins;

        return bx + by*_xbins + bz*_xbins*_ybins;
    }
}
