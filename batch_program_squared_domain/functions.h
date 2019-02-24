#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <string>

enum scalar_functions
{
    En,
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
enum weight
{
    area,
    inverse_distance_centroid,
    solid_angle
};
enum methods
{
    PCE,
    AGS,
    LSDD,
    LR,
    LR_centroids

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
    interior,
    boundaries
};


//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

std::string f_names[4]
{

    "axsin(3bx)sin(by^2)",
    "(b(x-1/2))^2+(c(y-1/2))^2",
    "asin(bx)cos(by)",
    "x+y"

};
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

std::string verts_names[3]
{

    "all",
    "interior",
    "baoundaries"

};
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

std::string weight_names[3]
{

    "Area",
    "InverseCentroidDistance",
    "SolidAngle"

};

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


std::string method_names[5]
{
    "PCE",
    "AGS",
    "LSDD",
    "LR",
    "LR_CENTROIDS"



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
    "Non Uniform",
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

std::string err_names[2]
{
  "Absolute",
  "Relative"
};


//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#endif // FUNCTIONS_H
