#include "lammps_data.h"
#include "string_tools.h"
#include "histogram.h"
#include <fstream>
#include <iostream>
#include <map>
#include <cmath>

using namespace std;

namespace cg {

    // Returns the scaled coordinate xs in [0,1).
    inline double to_scaled(double x, double dx) {
        x /= dx;
        while (x <  0.0) x += 1.0;
        while (x >= 1.0) x -= 1.0;
        return x;
    }

    // Computes the distance between two coordinates.
    inline double distance(const Coord &a, const Coord b) {    
        double dx=a.x-b.x, dy=a.y-b.y, dz=a.z-b.z;
        return sqrt(dx*dx+dy*dy+dz*dz);
    }

    // Reads the data from LAMMPS.
    LAMMPS_Data::LAMMPS_Data(string datapath, string dumppath) {
        if (datapath.size()) read_data(datapath);
        if (dumppath.size()) read_dump(dumppath);        
        type_counts();
    }

    // Counts how many of each atom type are present in the system.
    void LAMMPS_Data::type_counts() const {
        std::map<int, int> count;
        for (auto a=_system.atoms.cbegin(); a!=_system.atoms.cend(); ++a) {
            ++count[a->type];
        }        
        cout << "Type counts\n";
        for (auto t=count.cbegin(); t!=count.cend(); ++t)
            cout << t->first << ": " << t->second << "\n";        
    }

    // Reads a LAMMPS input file.
    void LAMMPS_Data::read_data(string path) {
        fstream fid(path);
        if (!fid) {
            cout << "Error opening output file " << path << "!\n";  
        }
        int section=0;    
        while (fid) {
            auto line = split(read_line(fid));
            if (line.empty()) continue;
            else if (line[0]=="Masses")    section=1;
            else if (line[0]=="Pair")      section=2;
            else if (line[0]=="Bond")      section=3;
            else if (line[0]=="Angle")     section=4;
            else if (line[0]=="Dihedral")  section=5;
            else if (line[0]=="Atoms")     section=6;
            else if (line[0]=="Bonds")     section=7;
            else if (line[0]=="Angles")    section=8;
            else if (line[0]=="Dihedrals") section=9;

            switch (section) {
            case 1: 
                if (line.size()==2) { 
                    //_system.atom_types.push_back(str2u32(line[0]));
                    //_system.masses.push_back(str2dbl(line[1]));
                }
                break;
            case 6:
                Atom newatom;
                if (line.size() >= 6) { 
                    newatom.molecule   = str2u32(line[1]);
                    newatom.type       = str2u32(line[2]);
                    _system.atoms.push_back(newatom);
                }

                break;
            case 7:
                if (_system.connect.empty()) {
                    _system.connect.assign(_system.atoms.size(), vector<int>(0));
                }
                if (line.size()==4) {
                    Bond newbond;
                    //newbond.id    = str2u32(line[0]);
                    //newbond.type  = str2u32(line[1]);
                    newbond.atom1 = str2u32(line[2])-1;
                    newbond.atom2 = str2u32(line[3])-1;
                    //_system.bonds.push_back(newbond);
                    _system.connect[newbond.atom1].push_back(newbond.atom2);
                    _system.connect[newbond.atom2].push_back(newbond.atom1);
                }
                break;
            case 8: break; // No need to keep going.
            default: continue;       
            };
        }
        // Report total number of bonds.
        int num_bonds = 0;
        for (auto c=_system.connect.cbegin(); c!=_system.connect.cend(); ++c) {
            num_bonds += c->size();
        }
        cout << "Read " << num_bonds/2 << " bonds from the input file.\n";            
    }

    // Reads a LAMMPS dump output file, in atom format.
    void LAMMPS_Data::read_dump(string path) {
        fstream fid(path);
        if (!fid) cout << "Error opening output file " << path << "!\n";
        else      cout << "Opened dump " << path << ".\n";

        SnapshotDump *step=NULL;
        while (fid) {
            auto line = trim(read_line(fid));
            if (line.empty()) continue;
            else if (line=="ITEM: TIMESTEP") {
                _dump.push_back(SnapshotDump());
                step = &_dump.back(); 
                step->timestep = str2u32(trim(read_line(fid)));
            }
            else if (line=="ITEM: NUMBER OF ATOMS") {
                if (!step) { 
                    cout << "Error: NUMBER OF ATOMS specified before TIMESTEP\n";
                    exit(1);
                }
                unsigned n = str2u32(trim(read_line(fid)));
                step->scaled_coordinates.assign(n, Coord());
            }
            else if (line.find("ITEM: BOX BOUNDS")==0) {
                if (!step) { 
                    cout << "Error: BOX BOUNDS specified before TIMESTEP\n";
                    exit(1);
                }
                auto xb  = split(read_line(fid));
                auto yb  = split(read_line(fid));
                auto zb  = split(read_line(fid));
                step->dx = fabs(str2dbl(xb.at(1))-str2dbl(xb.at(0)));
                step->dy = fabs(str2dbl(yb.at(1))-str2dbl(yb.at(0)));                                
                step->dz = fabs(str2dbl(zb.at(1))-str2dbl(zb.at(0)));
            }
            // The only thing left should be the ATOMS data.
            else {
                auto pieces = split(line);
                if (pieces.size() < 4) continue;
                if (pieces[0]=="ITEM:" && pieces[1]=="ATOMS") {
                    vector<string> var(pieces.begin()+2, pieces.end());
                    // Search for coordinate and tag columns.
                    step->xc = step->yc = step->zc = step->zc = -1;
                    step->xs = step->ys = step->zs = 0;
                    for (size_t i=0; i<var.size(); ++i) {
                        if (var[i]=="x")  {step->xc=i;}
                        if (var[i]=="xs") {step->xc=i; step->xs=1; }
                        if (var[i]=="y")  {step->yc=i;}
                        if (var[i]=="ys") {step->yc=i; step->ys=1; }
                        if (var[i]=="z")  {step->zc=i;}
                        if (var[i]=="zs") {step->zc=i; step->zs=1; }
                        if (var[i]=="id")  step->tc=i;
                    }
                    if (step->xc<0 || step->yc<0 || step->zc<0) {
                        cout << "Error: coordinate column not found\n";
                        exit(1);
                    }
                    if (step->tc<0) {
                        cout << "Error: atom tag column not found.\n";
                        exit(1);
                    }
                   continue;
                }
                if (!step) { 
                    cout << "Error: data encountered before TIMESTEP\n";
                    exit(1);
                }
                auto maxc=max(step->tc,max(step->xc,max(step->yc,step->zc)));
                auto minc=min(step->tc,min(step->xc,min(step->yc,step->zc)));
                if (maxc >= (int)pieces.size()) {
                    cout << "Not enough columns in dump data.\n";
                    exit(1);
                }
                if (minc < 0) {
                    cout << "Missing data column in dump data.\n";
                    exit(1);
                }
                unsigned id = str2u32(pieces[step->tc])-1;
                if (step->scaled_coordinates.size() <= id) {
                    cout << "Error: invalid atom id found " << id << "\n";
                    exit(1);
                }
                Coord &r = step->scaled_coordinates[id];
                r.x = str2dbl(pieces[step->xc]);
                r.y = str2dbl(pieces[step->yc]);
                r.z = str2dbl(pieces[step->zc]);

                if (!step->xs) r.x = to_scaled(r.x, step->dx);
                if (!step->ys) r.y = to_scaled(r.y, step->dy);
                if (!step->zs) r.z = to_scaled(r.z, step->dz);
            }
        }
    }

    // Returns the coordinate of an atom.
    double LAMMPS_Data::pair_distance(int step, int i, int j) const {
        const Coord &ri = _dump[step].scaled_coordinates[i];
        const Coord &rj = _dump[step].scaled_coordinates[j];

        double dx=fabs(ri.x-rj.x), dy=fabs(ri.y-rj.y), dz=fabs(ri.z-rj.z);  
        dx = min(dx, 1.0-dx)*_dump[step].dx;
        dy = min(dy, 1.0-dy)*_dump[step].dy;
        dz = min(dz, 1.0-dz)*_dump[step].dz;            
        return sqrt(dx*dx+dy*dy+dz*dz);       
    }

    // Returns pair distance between two particles (including periodic images).
    void LAMMPS_Data::pair_distances(int step, int i, int j, Histogram &h) const {
        auto dump       = _dump[step];
        const Coord &ri = dump.scaled_coordinates[i];
        const Coord &rj = dump.scaled_coordinates[j];
                
        double dx  = dump.dx*fabs(ri.x-rj.x);
        double dy  = dump.dy*fabs(ri.y-rj.y);
        double dz  = dump.dz*fabs(ri.z-rj.z);                
        double dx2[] = {dx, dx-dump.dx};  
        double dy2[] = {dy, dy-dump.dy};
        double dz2[] = {dz, dz-dump.dz};
        for (int a=0; a<2; ++a) {
            dx2[a] *= dx2[a];
            dy2[a] *= dy2[a];
            dz2[a] *= dz2[a];
        }
        
        for (int a=0; a<2; ++a) {
            for (int b=0; b<2; ++b)
                for (int c=0; c<2; ++c) 
                    h.add(sqrt(dx2[a] + dy2[b] + dz2[c]));
        }
    }

    // Returns the maximum possible distance between atoms in a periodic box.
    double LAMMPS_Data::min_box_size(int step) const {        
        const auto &d = _dump[step];
        return min(min(d.dx, d.dy), d.dz);       
    }

    //! Print atom position.
    void LAMMPS_Data::print_atom(int step, int i) const {
        const Coord &r = _dump[step].scaled_coordinates[i];
        cout << "["<<i<<"]: = ("<<r.x<<", "<<r.y<<", "<<r.z<<").\n";
    }

    //! Returns the volume at a time step.
    double LAMMPS_Data::volume(int step) const {
        auto s = _dump[step];
        return s.dx*s.dy*s.dz;
    }
}
