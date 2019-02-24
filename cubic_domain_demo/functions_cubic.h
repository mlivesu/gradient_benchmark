#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <string>

enum scalar_functions
{
    Lobb,
    Paraboloid,
    Eggs,
    Test,
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
    LRCentroids,
    LSDDCentroids,
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
    non_uniform,
    anisotropic,
    grid
};

enum type_of_vertices
{
    all,
    no_srf,
    only_srf

};

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

std::string f_names[4]
{

    "Lobb Function",
    "Paraboloid",
    "Eggs",
    "Test"

};

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


std::string method_names[6]
{
    "PCE",
    "AGS",
     "LSDD",
    "LR",
    "LR_Centroids",
    "LSDD_Centroids"



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

std::string tri_names[4]
{
    "Regular",
    "Delaunay",
    "Non-Uniform",
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
