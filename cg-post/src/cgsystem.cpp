#include "cgsystem.h"
#include "lammps_data.h"
#include <algorithm>
#include <numeric>
#include <cmath>
#include "string_tools.h"
#include "bintable.h"

using std::count;

namespace cg { 
    // Constructor for the cg_system.
    CgSystem::CgSystem(string in) 
        : _lammps_data(NULL) {

        // Set default histogram sizes.
        _adf_param[0] = 0.0;       
        _adf_param[1] = 180.0;
        _adf_param[2] = 1.0;

        _rdf_param[0] = 0.0;
        _rdf_param[1] = 20.0;
        _rdf_param[2] = 0.2;
        
        _bdf_param[0] = 0.0;
        _bdf_param[1] = 15.0;
        _bdf_param[2] = 0.1;

        read_input(in);

        for (auto &bead_iter: _bead_definitions) {
            unsigned old_size = _beads.size();
            find_beads(bead_iter.first);
            cout << "Bead " << bead_iter.second.atoms
			 	<< "\t" << _beads.size() - old_size << " found.\n";
        }
        find_bead_neighbors();
    }
    // Clears the LAMMPS data.
    CgSystem::~CgSystem() { delete _lammps_data; }    

    // Returns the number of timesteps.
    int CgSystem::num_timesteps() const { return _lammps_data->_dump.size(); }

    // Returns the element of an atom type.
	string CgSystem::type_to_element(int type) const { 
        return _type.find(type)->second; 
    }

    // Reads the input file for the cg-analysis.
    void CgSystem::read_input(string path)  {
        string lmp_data, lmp_dump;
        std::fstream fid(path.c_str(), std::ios::in);
        if (!fid) {
            cout << "Error input " << path << " cannot be opened.\n";
            exit(1);
        }
        while (fid) {
            auto line = read_line(fid);
            if (line.empty()) continue;
            line = line.substr(0, line.find_first_of("#"));
            auto cmd = split(line);
            const size_t n = cmd.size();

            if (n==0) /* empty line */;
            else if (cmd[0]=="input"  && n==2) lmp_data = cmd[1];
            else if (cmd[0]=="dump"   && n==2) lmp_dump = cmd[1];        
            else if (cmd[0]=="type"   && n==3) _type[str2u32(cmd[1])] = cmd[2];
            else if (cmd[0]=="output" && n==2) _output_tag = cmd[1];
            else if (cmd[0]=="rdf" || cmd[0]=="bdf"||cmd[0]=="adf") {
                if (n != 4) {
                    cout << "Error: " << cmd[0] << " beg end inc " << "\n.";
                    exit(1);
                }
                
                double *p;
                if      (cmd[0]=="adf")  p = _adf_param;
                else if (cmd[0]=="rdf")  p = _rdf_param;
                else p = _bdf_param;                
                
                p[0] = str2dbl(cmd[1]); 
                p[1] = str2dbl(cmd[2]);
                p[2] = str2dbl(cmd[3]);
            }            
            else if (cmd[0]=="bead"   && n>2) {
                auto &bead_id = cmd[1];
                _bead_definitions[bead_id] = BeadDefinition();
                for (auto e=cmd.cbegin()+2; e!=cmd.cend(); ++e) {
                    _bead_definitions[bead_id].atoms.push_back(*e);
                }
            }
            else if (cmd[0]=="share" && n>2) {
                auto &bead_id  = cmd[1];
                auto bead_iter = _bead_definitions.find(bead_id);

                if (bead_iter == _bead_definitions.end()) {
                    cout << "Error: bead " << bead_id << " not yet defined.\n";
                    exit(1);
                }
                // Push back any atoms indices after the id to indicate sharing.
                for (auto i=cmd.cbegin()+2; i!=cmd.cend(); ++i)
                   bead_iter->second.share.push_back(str2u32(*i));
            }
            // Must be called after a bead command.
            else if (cmd[0]=="center" && n==3) {
                auto &bead_id  = cmd[1];
                auto bead_iter = _bead_definitions.find(bead_id);
                if (bead_iter == _bead_definitions.end()) {
                    cout << "Error: bead " << bead_id << " not yet defined.\n";
                    exit(1);
                }
                bead_iter->second.center = str2u32(cmd[2]);
            }
            else cout << "Invalid command \"" << join(cmd) << "\".\n";
        }
        
        _lammps_data = new cg::LAMMPS_Data(lmp_data, lmp_dump);
        _used_atoms.assign(_lammps_data->_system.atoms.size(), 0);
    }

    // Crawls though the connectivity and determines beads.
    void CgSystem::find_beads(const std::string &beadtype)
    {
        cout << "Finding beads of type " << beadtype << "\n";
		auto &beadDef = _bead_definitions[beadtype];
        auto &share   = beadDef.share;
        const AtomicSystem &sys = _lammps_data->_system;
        if (beadDef.atoms.empty()) {
            cout << "CG bead definition has no atoms.\n";
            return;
        }       
        string start_type = beadDef.atoms.front();
        // Loop over atoms try to find first bead atom.
        for (size_t i=0; i<sys.atoms.size(); ++i) {
            // Skip if atom is already in a bead and not shared.
            bool shared = count(share.cbegin(), share.cend(), 0)>0;                                        
            if (_used_atoms[i] && !shared) continue;
            if (!match(type_to_element(sys.atoms[i].type), start_type)) continue;

            // Start potential bead at i.        
            std::vector<int> pot_bead(1, i);        
            while (true) {
                auto length = pot_bead.size();
                // Bead starts at i.
                if (length == beadDef.atoms.size()) {
                    Bead new_bead;
                    new_bead.type = beadtype;
                    new_bead.center = pot_bead[beadDef.center];
                    new_bead.ends.push_back(pot_bead.front());
                    // If bead is length = 1, don't repeat end atom.
                    if (length > 1) new_bead.ends.push_back(pot_bead.back());                    
                    _beads.push_back(new_bead);
                    for (auto id=pot_bead.cbegin(); id!=pot_bead.cend(); ++id) 
                        ++_used_atoms[*id];
                    break;
                }

                int current      = pot_bead.back();
                string next_type = beadDef.atoms[length];
                shared = count(share.cbegin(),share.cend(),length)>0;
                                        
                // Try to find next bead.
                auto next = sys.connect[current].cbegin();            
                for (; next != sys.connect[current].cend(); ++next) {
                    if (_used_atoms[*next] && !shared) continue;
					// skip if focus atom already in potential bead
                    if (count(pot_bead.begin(), pot_bead.end(), *next)) continue;
                    if (match(type_to_element(sys.atoms[*next].type), next_type)) {
                        pot_bead.push_back(*next);
                        break;
                    }
                }
                // No bead starts at i.
                if (next == sys.connect[current].cend()) break;
            }        
        }        
    }

    // Finds the neighbors for each bead.
    void CgSystem::find_bead_neighbors()  {        
        const AtomicSystem &sys = _lammps_data->_system;                 
        int ct=0;
        for (auto i=_beads.begin(); i!=_beads.end()-1; ++i) {
            // Check if bead j is bonded to bead i.
            for (auto j=i+1; j!=_beads.end(); ++j) {
                // Check if any end atoms on i are bonded to any end atoms on j
                // or if any end atoms of i are also end atoms of j.
                for (auto ie=i->ends.cbegin(); ie!=i->ends.cend(); ++ie) {
                    if (count(j->ends.cbegin(), j->ends.cend(), *ie)) {
                        i->neighbors.push_back(j);
                        j->neighbors.push_back(i);
                        ++ct;
                        break;
                    }
                    const auto &iconn = sys.connect[*ie];
                    for (auto je=j->ends.cbegin(); je!=j->ends.cend(); ++je) {
                        if (count(iconn.cbegin(), iconn.cend(), *je)) {
                            i->neighbors.push_back(j);
                            j->neighbors.push_back(i);
                            ++ct;
                            break;
                        }
                    }
                    if (count(i->neighbors.cbegin(),i->neighbors.cend(),j)) 
                        break;
                }
            }
        }        
        cout << "Found " << ct << " number of bead neighbors.\n";
    }

    // Computes the bond-length distrubution function between beads.
    Histogram CgSystem::compute_bdf(string b1, string b2, int step) const {
        Histogram bdf(_bdf_param[0], _bdf_param[1], _bdf_param[2]);        
        
        if (b1 > b2) std::swap(b1, b2);
        // Loop over all pairs of beads (w/ no cutoff).
        for (auto i=_beads.cbegin(); i!=_beads.cend(); ++i) {            

            // If bead i matches b1 or b2, then bead j must match the other.
            string bj;
            if      (i->type == b1) bj = b2;
            else if (i->type == b2) bj = b1;
            else continue;

            for (auto j=i+1; j!=_beads.cend(); ++j) {                
                if (!neighbors(i, j)) continue;
                if (j->type != bj) continue;

                double d = _lammps_data->pair_distance(step, i->center, j->center);
                bdf.add(d);
            }        
        }
        return bdf;
    }

    // Computes the bond-angle distrubution function between beads.
    Histogram CgSystem::compute_adf(string b1, string b2, string b3, int step) const {
        static double rad2deg = 180.0/acos(-1.0);
        Histogram adf(_adf_param[0], _adf_param[1], _adf_param[2]);            
        // Loop over all pairs of beads (w/ no cutoff).
        for (auto i=_beads.cbegin(); i!=_beads.cend(); ++i) {
            if (i->neighbors.size() < 2) continue;
            int j=i->neighbors[0]->center, k=i->neighbors[1]->center;
            // For now we'll just assume always 2 neighbors.
            double a = _lammps_data->pair_distance(step, i->center, j);
            double b = _lammps_data->pair_distance(step, i->center, k);
            double c = _lammps_data->pair_distance(step, j, k);
            double q = acos((a*a + b*b - c*c)/(2.0*a*b));
            adf.add(q*rad2deg);
        }
        return adf;
    }

    // Computes the radial distribution function between beads.
    Histogram CgSystem::compute_rdf(string b1, string b2, int step) const {

        Histogram rdf(_rdf_param[0], _rdf_param[1], _rdf_param[2]);

        // Loop over all pairs of beads (w/ no cutoff).
            // If bead i matches b1 or b2, then bead j must match the other.
        for (auto i=_beads.cbegin(); i!=_beads.cend(); ++i) {
            string bj;
            if      (i->type == b1) bj = b2;
            else if (i->type == b2) bj = b1;
            else continue;
                
            for (auto j=i+1; j!=_beads.cend(); ++j) {
                if (j->type != bj || neighbors(i,j)) continue;
                _lammps_data->pair_distances(step, i->center, j->center, rdf);            
            }
        }

        if (b1 == b2) {
            double n   = (double)bead_count(b1);
            double rho = number_density(b1, step);
            double wt  = 2.0/(n*rho);
            rdf.scale(wt, true);
        }
        else {
            double n1   = (double)bead_count(b1);
            double n2   = (double)bead_count(b2);
            double rho1 = number_density(b1, step);
            double rho2 = number_density(b2, step);
            double wt   = 2.0 / (n1*rho2 + n2*rho1);
            rdf.scale(wt, true);
        }                        
        return rdf;
    }

    // Determines if beads, i, j are 1st neighbors.
    bool CgSystem::neighbors(BeadCIter i, BeadCIter j) const {
        return std::count(i->neighbors.cbegin(), i->neighbors.cend(), j)>0;
    }

    bool CgSystem::neighbors(int i, int j) const {
        return std::count(_beads[i].neighbors.cbegin(), 
                          _beads[i].neighbors.cend(), _beads.cbegin() + j);
    }

    // Determines if beads, i, j share a common neighbor.
    bool CgSystem::second_neighbors(BeadCIter i, BeadCIter j) const {
        for (auto k=i->neighbors.cbegin(); k!=i->neighbors.cend(); ++k) {
            if (std::count(j->neighbors.cbegin(), j->neighbors.cend(), *k)) return true;
        }
        return false;
    }

    // Returns the number of beads of a given type.
    int CgSystem::bead_count(string bead_type) const {
        int count = 0;
        for (auto b=_beads.begin(); b!=_beads.end(); ++b) {
            count += int(b->type == bead_type);
        }        
        return count;
    }

    // Returns the number of beads divided by the system volume.
    double CgSystem::number_density(string bead_type, int step) const {
        double count = double(bead_count(bead_type));
        double vol   = _lammps_data->volume(step);
        return count/vol;        
    }

    // Returns the x y and z coordinates of a bead.
    Coord CgSystem::bead_position(int step, int i) const {
        // For now we just use center atom.
        int c = _beads[i].center;
        return _lammps_data->atom_position(step, c);
        
    }

    //! Returns a vector of all bead types defined.
    std::vector<string> CgSystem::defined_bead_types() const {
        std::vector<string> types;
        for (auto &it: _bead_definitions) types.push_back(it.first);
        return types;
    }

    // Checks if type1 matches type2 (allows wildcards on type2).
    bool match(std::string type1, std::string type2) {
        if (type1.find('*') != type1.npos) {
            std::cout << "Warning: Wildcard found in type.\n";
        }        
        
        if (type2.find('*') == type2.npos) return type1 == type2;
            
        // If there is a wildcard, then we need to match everything before the wildcard.
        auto pieces = split(type2, "*");
        // Type2 was only stars.
        if (pieces.empty()) return true; 

        size_t loc=type1.find(pieces[0]);
        if (loc > 0) return false;
        for (auto p=pieces.cbegin()+1; p!=pieces.cend(); ++p) {
            loc = type1.find(*p, loc+1);
            if (loc == type1.npos) return false;
        }
        return true;
    }

    double CgSystem::box_dx(int step) const {
        return _lammps_data->box_dx(step);
    }
    double CgSystem::box_dy(int step) const {
        return _lammps_data->box_dy(step);
    }
    double CgSystem::box_dz(int step) const {
        return _lammps_data->box_dz(step);
    }
}

