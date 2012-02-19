#include "lammps_data.h"
#include "string_tools.h"
#include "cgsystem.h"
#include <string>
#include <vector>

using namespace std;

int main(int n, char **argv) 
{    
    string input = "cg.ini";
    vector<string> args(argv, argv+n);
    for (auto a=args.cbegin()+1; a!=args.cend(); ++a) {
        if (*a=="-i"&&++a!=args.cend()) input = *a;
        if (*a=="--test")               string_test();
    }    
    cg::CgSystem model(input);

    auto bead_types = model.defined_bead_types();
    
    cout << "Computing histograms for " << model.num_timesteps() << " steps.\n";
    for (auto b1=bead_types.cbegin(); b1!=bead_types.cend(); ++b1) {
        for (auto b2=b1; b2!=bead_types.cend(); ++b2) {
            cg::Histogram rdf, bdf, adf;
            cout << "Computing " << *b1 << ":" << *b2 << " correlations.\n";
            #pragma omp parallel for
            for (int ts=0; ts<model.num_timesteps(); ++ts) {                    
                auto new_rdf = model.compute_rdf(*b1, *b2, ts);
                auto new_bdf = model.compute_bdf(*b1, *b2, ts);

                //adf += model.compute_adf(*b1, *b2, *b3, ts);
                #pragma omp critical 
                rdf += new_rdf;
                #pragma omp critical 
                bdf += new_bdf;
            }

            if (rdf.histogram.size()==0) return 0;
            rdf.scale(1.0/model.num_timesteps(), false);    
            bdf.scale(1.0/model.num_timesteps(), false);
            //adf.scale(1.0/model.num_timesteps(), false);

            bdf.normalize();
            //adf.normalize();

            string tag   = model.output_tag();
            rdf.write("rdf-"+tag+"_"+*b1+*b2+".txt");
            bdf.write("bdf-"+tag+"_"+*b1+*b2+".txt");
            //adf.write("adf-"+tag+*b1+*b2+*b3+".txt");
        }        
    }
}
