#pragma once
#include <vector>
#include <string>
#include <map>
#include "histogram.h"
#include "coord.h"

namespace cg {
using std::string;

class LAMMPS_Data;

//! Define beads as a sequence of atom identifiers (string)
struct BeadDefinition {
    //! Atoms that make up the bead.
    std::vector<string> atoms;
    //! Index of atom specifying the bead center.
    int center;
    //! Atoms that can be shared by two beads.
    std::vector<int> share;           
};

struct Bead;
typedef std::vector<Bead> VBead;
typedef VBead::const_iterator BeadCIter;

//! Object that defines a coarse grain bead.
struct Bead {
    //! Type of this bead.
    string type;                 
    //! Center atom of bead.
    int center;
    //! End atoms of bead.
    std::vector<int> ends;
    //! Neighbors to this bead (probably 2).    
    std::vector<BeadCIter> neighbors;
};

class CgSystem {
public:
    //! Initializes the cg_system.
    CgSystem(string in);
    ~CgSystem();

    //! Reads the input file to build the cg system.
    void read_input(string path);
    //! Returns the number of timesteps.
    int num_timesteps() const;
    //! Returns the element of an atom type.
	string type_to_element(int type) const;
    //! Finds the atoms that make up the bead centers.
    void find_beads(const string &beadtype);    
    //! Finds the neighbors for each bead.
    void find_bead_neighbors();
    //! Computes the bond-length distrubution function between beads.
    Histogram compute_bdf(string b1, string b2, int step) const;
    //! Computes the bond-angle distrubution function between beads.
    Histogram compute_adf(string b1, string b2, string b3, int step) const;
    //! Computes the radial distribution function between beads.
    Histogram compute_rdf(string b1, string b2, int step) const;
    //! Determines if beads, i, j are 1st neighbors.
    bool neighbors(BeadCIter i, BeadCIter j) const;
    bool neighbors(int i, int j) const;
    //! Determines if beads, i, j are 2nd neighbors.
    bool second_neighbors(BeadCIter i, BeadCIter j) const;
    //! Returns the output tag.
    string output_tag() const { return _output_tag; }
    //! Returns the number of beads by type.
    int bead_count(string bead_type) const;
    //! Returns the average number density of the system.
    double number_density(string bead_type, int step) const;
    //! Returns the number of beads defined in the system.
    int num_bead_types() const { return _bead_definitions.size(); }
    //! Returns a vector of all bead types defined.
    std::vector<string> defined_bead_types() const;
    //! Returns the x y and z coordinates of a bead.
    Coord bead_position(int step, int i) const;
    int num_beads() const { return _beads.size(); }

    double box_dx(int step) const;
    double box_dy(int step) const;
    double box_dz(int step) const;
    
private:
    //! Definition of each bead type.
    std::map<string, BeadDefinition> _bead_definitions;

    //! Atoms not yet mapped to beads.
    std::vector<char> _used_atoms;
    //! Beads defined in this system.
    std::vector<Bead> _beads;
    //! Tag to add to the output files.
    string _output_tag;
    //! Converts from integer type (LAMMPS) to element.
    std::map<int, string> _type;
    LAMMPS_Data *_lammps_data;

    //! Historgram parameters are {start, end, increment}.
    double _rdf_param[3], _bdf_param[3], _adf_param[3];
};
    //! Checks if type1 matches type2 (allows wildcards on type2).
    bool match(string type1, string type2);

}
