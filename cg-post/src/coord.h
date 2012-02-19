#pragma once

//! Position in three dimensions. 
struct Coord {
    Coord() {}
    Coord(double x1, double x2, double x3) : x(x1), y(x2), z(x3) {}
    double x,y,z;
};

