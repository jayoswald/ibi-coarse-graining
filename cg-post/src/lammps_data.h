#pragma once
#include <string>
#include <vector>

namespace cg {
using std::string;

//! Position of an atom.
struct Coord {
    Coord() {}
    Coord(double x1, double x2, double x3) : x(x1), y(x2), z(x3) {}
    double x,y,z;
};

//! Output from LAMMPS dump style atom.
struct SnapshotDump {
    //! Actual count of timestep.
    int timestep; 
    //! Box sizes along each dimension.
    double dx, dy, dz;  
    //! Indices where the required data is located.
    int xc, yc, zc, tc; 
    //! Whether each dof is scaled or not.
    bool xs, ys, zs;
    std::vector<Coord> scaled_coordinates;
};

class Histogram;

//!
struct Atom {
    int type;           // Atom type - one-based type index.
    int molecule;       // Molecule atom is part of - one based.
};

//! Information stored for each bond between two atoms.
struct Bond {
    int id, type, atom1, atom2;
};

//! 
struct AtomicSystem {    
    //std::vector<double> masses;           //! Mass of each atom. (not needed)
    std::vector<Atom> atoms;                //! Each atom stored at index id-1.
    //std::vector<bond> bonds;
    std::vector<std::vector<int>> connect;  //! Atom connectivity - links ids.
};

class LAMMPS_Data {
    friend class CgSystem;
public:
    //! Builds lammps data structure; optional: read input and dump.
    LAMMPS_Data(string datapath="", string dumppath="");
    //! Reports how many of each atom type are present in the system.
    void type_counts() const;
    //! Reads a LAMMPS input file.
    void read_data(string path);
    //! Reads a LAMMPS dump output file.
    void read_dump(string path);
    //! Returns the coordinate of an atom.
    double pair_distance(int step, int i, int j) const;
    //! Returns pair distance between two particles (including periodic images).
    void pair_distances(int step, int i, int j, Histogram &h) const;
    //! Returns the maximum possible distance between atoms in a periodic box.
    double min_box_size(int step) const;
    //! Print atom position.
    void print_atom(int step, int i) const;

    //! Returns the volume at a time step.
    double volume(int step) const;
private:
    AtomicSystem _system;    
    std::vector<SnapshotDump> _dump;
};
}

