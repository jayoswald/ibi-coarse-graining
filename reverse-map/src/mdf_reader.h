#pragma once
#include <string>
#include <map>
#include <vector>

typedef std::map<std::string, std::vector<int>> 


class Atom {
public:

private:
    std::string _name, _element, _fftype;
    double _charge;
};

class InsightData {
public:
    // Reads Insight .car and .mdf files.
    class InsightData(std::string basename);
    
private:

    std::string mdf_header;
    std::string car_header;
};


