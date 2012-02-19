#pragma once
#include <vector>

namespace cg {

    class CgSystem;
    //! Constructs a bin table for fast neighbor searches.
    class BinTable {
    public:
        //! Constructs a bin table from an output dump.
        BinTable(double cutoff, int step, const CgSystem *cg); 
        
        //! Returns an array of all atoms that can be within cutoff of i.
        std::vector<const std::vector<int>*> neighbors(int i) const;
        
    private:
        //! Returns the bin number of atom i.
        int _find_bin(int i) const;

        // Sets the number of bins per cutoff.
        int _granularity;
        int _xbins, _ybins, _zbins;
        std::vector<std::vector<int>> _bins;

        //! Timestep to pull data from.        
        int _step;
        //! Cutoff radius.
        double _cutoff;
        //! System to get coordinates from.
        const CgSystem *_cgsys;
    };
}

