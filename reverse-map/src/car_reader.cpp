#include "car_reader.h"
#include <fstream>
#include <string>
#include <iostream>
#include "string_tools.h"

using namespace std;

void car_reader(const string &path) {
    fstream fid(path, ios::in);
    if (!fid.is_open()) cout << "Cannot open car file: " << path << "\n"; 

    bool pbc = false;
    while (fid) { 
        auto line   = trim(read_line(fid));
        auto pieces = split(line);
        if (pieces.empty() || line == "end") {
            // Empty line.
        }
        else if (pieces[0][0] == '!') {
            // Do something with comment string.
        } 
        else if (startswith(line, "Materials Studio Generated")) {
            // Also a comment, ignore.
        }
        else if (startswith(line, "PBC=")) {
            pbc = line == "PBC=ON";
        }
        else if (pieces.size() == 9) {
            string atom=pieces[0], molecule=pieces[4], fftype=pieces[6];
            string element=pieces[7];
            double x = from_string<double>(pieces[1]);
            double y = from_string<double>(pieces[2]);
            double z = from_string<double>(pieces[3]);
            double q = from_string<double>(pieces[8]);
        }
        // Don't know what to do.
        else cout << line << "\n";
    }
}

