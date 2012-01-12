#pragma once
#include <vector>

namespace cg {

    class SnapshotDump;
    //! Constructs a bin table for fast neighbor searches.
    class BinTable {
        //! Constructs a bin table from an output dump.
        BinTable(double cutoff, const SnapshotDump &dump);
        
        //! Returns an array of all atoms that can be within cutoff of i.
        std::vector<int> neighbors(int i) const;
        
    private:
        // Sets the number of bins per cutoff.
        static const int _granularity = 2;
        //! Returns the bin number of atom i.
        int _find_bin(int i) const;

        int _xbins, _ybins, _zbins;
        std::vector<std::vector<int>> _bins;
        const SnapshotDump* _dump;

        double _cutoff;
    };


}

