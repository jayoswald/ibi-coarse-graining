#include "bintable.h"
#include "lammps_data.h"

namespace cg {

    BinTable::BinTable(double cutoff, const SnapshotDump &dump) {
        // Keep pointer reference to dump so we can look up bins later.
        _dump = &dump;
        _cutoff = cutoff;
        // Number of bins along each direction.
        // Bins are sized such that two atoms within 2 neighbor bins
        // must be within the cutoff distance.
        _xbins = int(_granularity*dump.dx/cutoff);
        _ybins = int(_granularity*dump.dy/cutoff);
        _zbins = int(_granularity*dump.dy/cutoff);
        
        _bins.assign(_xbins*_ybins*_zbins, std::vector<int>());

        for (int i=0; i<dump.scaled_coordinates.size(); ++i) {
            _bins[_find_bin(i)].push_back(i);
        }
    }

    // Returns an array of all atoms that can be within cutoff of i.
    // This might not be the best implimentation because a lot of 
    // vectors need to be copied.
    std::vector<int> BinTable::neighbors(int i) const {
        // Scaled coordinate of the atom.
        auto xi = _dump->scaled_coordinates[i];

        int bi  = _find_bin(i);
        std::vector<int> neigh;

        int bx = bi % _xbins;
        int by = (bi/_xbins)%_ybins;
        int bz = (bi/_xbins/_ybins);


        int N = _granularity;
        for (int i=-N; i<=N; ++i) {
            int nx = bx + i;
            if (nx<0 || nx>=_xbins) continue;
            for (int j=-N; j<=N; ++j) {                
                int ny = by + j;
                if (ny<0 || ny>=_ybins) continue;
                for (int k=-N; k<=N; ++k) {
                    int nz = bz + k;
                    if (nz<0 || nz>=_zbins) continue;

                    // See if bin ijk should be added.
                    double dx = double(i-1)/double(_xbins) * _dump->dx;
                    double dy = double(j-1)/double(_ybins) * _dump->dy;
                    double dz = double(k-1)/double(_zbins) * _dump->dz;
                    if (dx*dx + dy*dy + dz*dz > _cutoff*_cutoff) continue;

                    int nb = nx + ny*_xbins + nz*_xbins*_ybins;

                    neigh.insert(neigh.end(), _bins[nb].begin(), _bins[nb].end());
                }
            }
        }
        return neigh; 
    }

    int BinTable::_find_bin(int i) const {
        auto xs = _dump->scaled_coordinates[i];
        int bx  = double(_xbins)*xs.x;
        int by  = double(_ybins)*xs.y;
        int bz  = double(_zbins)*xs.z;

        // Make sure we don't overflow if xs.? was 1.0.
        if (bx == _xbins) bx--;
        if (by == _ybins) by--;
        if (bz == _zbins) bz--;

        return bx + by*_xbins + bz*_xbins*_ybins;
    }
}
