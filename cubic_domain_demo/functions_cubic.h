#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <string>

enum scalar_functions
{
    FUNC_0,
    FUNC_1,
    FUNC_2,
    N_FUNCTIONS
};
enum error_type
{
    absolute,
    relative
};
enum mode
{
    total,
    norm,
    angle
};

enum methods
{
    PCE,
    AGS,
    LSDD,
    LR,
    N_METHODS

};

enum canvas
{
    I,
    II,
    III,
    IV
};
enum tri
{
    regular,
    Delaunay,
    anisotropic,
    N_tri
};

enum type_of_vertices
{
    all,
    no_srf,
    only_srf

};

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

std::string f_names[3]
{

    "Eggs",
    "Paraboloid",
    "Test"

};

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


std::string method_names[4]
{
    "PCE",
    "AGS",
     "LSDD",
    "LR"


};
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

std::string canvas_names[4]
{
    "I",
    "II",
    "III",
    "IV"

};
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

std::string tri_names[3]
{
    "Regular",
    "Delaunay",
    "Anisotropic"

};
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

std::string mode_names[3]
{
    "Total",
    "Norm",
    "Angle"


};
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

std::string vertices_names[3]
{
    "All vertices",
    "Exclude Srf",
    "Only Srf"


};
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

std::string err_names[2]
{
  "Absolute",
  "Relative"
};


//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#endif // FUNCTIONS_H
