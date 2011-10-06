#pragma once
#include <vector>
#include <string>
#include <map>
#include "histogram.h"

namespace cg {

class LAMMPS_Data;

//! Define beads as a sequence of atom identifiers (string)
struct BeadDefinition {
    std::vector<std::string> atoms;
    int center;
    //! Atoms that can be shared by two beads.
    std::vector<int> share;           
};

struct Bead;
typedef std::vector<Bead> VBead;
typedef VBead::const_iterator BeadCIter;

//! Object that defines a coarse grain bead.
struct Bead {
    int type;                         //! Type of this bead.
    int center;                       //! Center atom of bead.
    std::vector<int> ends;            //! End atoms of bead.
    std::vector<BeadCIter> neighbors; //! Neighbors to this bead (probably 2).    
};

class CgSystem {
public:
    //! Initializes the cg_system.
    CgSystem(std::string in);
    ~CgSystem();

    //! Reads the input file to build the cg system.
    void read_input(std::string path);
    //! Returns the number of timesteps.
    int num_timesteps() const;
    //! Returns the element of an atom type.
	std::string type_to_element(int type) const { return _type.find(type)->second; }
    //! Finds the atoms that make up the bead centers.
    void find_beads(size_t beadtype);    
    //! Finds the neighbors for each bead.
    void find_bead_neighbors();
    //! Computes the bond-length distrubution function between beads.
    Histogram compute_bdf(int b1, int b2, int step) const;
    //! Computes the bond-angle distrubution function between beads.
    Histogram compute_adf(int b1, int b2, int b3, int step) const;
    //! Computes the radial distribution function between beads.
    Histogram compute_rdf(int b1, int b2, int step) const;
    //! Determines if beads, i, j are 1st neighbors.
    bool neighbors(BeadCIter i, BeadCIter j) const;
    //! Determines if beads, i, j are 2nd neighbors.
    bool second_neighbors(BeadCIter i, BeadCIter j) const;
    //! Returns the output tag.
    std::string output_tag() const { return _output_tag; }
    //! Returns the number of beads by type.
    int bead_count(int type) const;
    //! Returns the average number density of the system.
    double number_density(int type, int step) const;

    //! Returns the number of beads defined in the system.
    int num_bead_types() const { return _bead_defs.size(); }
    
private:
    //! Definition of each bead type.
    std::vector<BeadDefinition> _bead_defs;
    //! Atoms not yet mapped to beads.
    std::vector<char> _used_atoms;
    //! Beads defined in this system.
    std::vector<Bead> _beads;
    //! Tag to add to the output files.
    std::string _output_tag;
    //! Converts from integer type (LAMMPS) to element.
    std::map<int, std::string> _type;
    LAMMPS_Data *_lammps_data;

    //! Historgram parameters are {start, end, increment}.
    double _rdf_param[3], _bdf_param[3], _adf_param[3];
};
    //! Checks if type1 matches type2 (allows wildcards on type2).
    bool match(std::string type1, std::string type2);

}
