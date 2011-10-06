#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>

using std::string;
using std::vector;
using std::fstream;
using std::cout;


//! Writes a vector to the console.
template<class T>
std::ostream& operator<<(std::ostream& o, std::vector<T> x)
{
    o << "{";
    if (x.size()) o << x.front();
    for (size_t i=1; i<x.size(); ++i)  o << ", " << x[i];
    return o << "}\n";
}

//! Reads the next line of a file.
inline std::string read_line(std::fstream &fid)
{
    std::string line;
    getline(fid, line);
    return line;
}

//! Converts a string to anything.
template<class T>
inline T from_string(const std::string &s, int *fail=NULL)
{
    std::istringstream ss(s);
    T t;
    ss >> t;
    if (ss.fail()) {        
        if (fail) *fail=1;
        else std::cout<<"Warning: bad string conversion:" + s + "\n";
        return T();
    }
    if (fail) *fail=0;
    return t;
}
//! Shortcut to convert strings to unsigned.
inline unsigned str2u32(const std::string &s, int *fail=0)
{
    return from_string<unsigned>(s,fail);
}
//! Shortcut to convert strings to double.
inline double str2dbl(const std::string &s, int *fail=0)
{
    return from_string<double>(s,fail);
}

//! Returns a copy of the string without trailing/preceding whitespace.
inline std::string trim(std::string s, std::string ws=" \t\n\r")
{
    if (s.empty()) return s;
    size_t a=s.find_first_not_of(ws);
    size_t b = s.find_last_not_of(ws)+1;
    return s.substr(a, b-a);
}

//! Splits a string like the python function.
inline std::vector<std::string> split(std::string s, std::string delims=" \t\n")
{
    std::vector<std::string> pieces;
    size_t begin=0, end=0;
    while (end != std::string::npos) {
        begin = s.find_first_not_of(delims, end);
        end   = s.find_first_of(delims, begin);
        if (begin != std::string::npos)
            pieces.push_back(s.substr(begin, end-begin));
    }
    return pieces;
}

//! Joins a vector, string back together.
inline std::string join(const std::vector<std::string> &v)
{
    if (v.empty()) return "";
    std::string s(v.front());
    std::vector<std::string>::const_iterator viter=v.begin()+1;
    for (; viter != v.end(); ++viter) s += " " + *viter;
    return s;
}

///////////////////////////////////////////////////////////////////////////////
// Converts anything to a string
///////////////////////////////////////////////////////////////////////////////
template<typename T>
inline string to_string(const T& x, int width=-1, int precision=-1) {
    std::ostringstream o;
    o<<x;
    return o.str();
}

//! Temporary main for testing/demoing the various functions.
static int string_test() {
    const char *test = "\tThis is a sentence with seven words.\n ";
    std::cout <<'"'<< test        << "\"\n";
    std::cout <<'"'<< trim(test)  << "\"\n";
    std::cout <<'"'<< split(test) << "\"\n";
    std::cout <<'"'<< join(split(test)) << "\"\n";
    std::cout << "42 is " << str2u32("forty-two") <<"\n";

    int fail;
    std::cout << "42 is " << str2u32("forty-two",&fail) <<"\n";
    if (fail) std::cout <<"Detected a bad conversion\n";

    std::vector<std::string> parts = split(test);
    for (size_t i=0; i<parts.size(); ++i)
        if (parts[i] == "sentence") std::cout << "found it\n";
    return 0;
}

