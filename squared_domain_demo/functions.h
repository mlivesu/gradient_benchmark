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
    LR_CENTROIDS

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
    anisotropic
};

enum vertices
{
    all,
    exclude_boundary,
    only_boundary
};

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

std::string f_names[3]
{

    "asin(bx)cos(by)",
    "b(x-1/2)^2+c(y-1/2)^2",
    "x+y"

};

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


std::string method_names[5]
{
    "PCE",
    "AGS",
     "LSDD",
    " LR",
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
    "All",
    "Exlude Boundary",
    "Only "


};
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

std::string err_names[2]
{
  "Absolute",
  "Relative"
};
double eval_func(const double x, const double y, const int func)
{
    switch(func)
    {
        case FUNC_0 : return x*sin(x)*sin(pow(1.5,y));
        case FUNC_1 : return x-y;
        default: assert(false && "Unknown Function");
    }
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#endif // FUNCTIONS_H
