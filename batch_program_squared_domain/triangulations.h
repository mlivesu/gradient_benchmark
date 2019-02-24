#ifndef TRIANGULATIONS_H
#define TRIANGULATIONS_H

#include <cinolib/triangle_wrap.h>
#include <cinolib/meshes/meshes.h>
#include<cinolib/drawable_vector_field.h>
#include <math.h>

using namespace cinolib;

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void anisotropic_tri(DrawableTrimesh<> & m, const double max_area=0.01, const double min_angle=10, const double y_anisotropy=5);
void regular_tri (DrawableTrimesh<> & m, int N);
void constraint_tri(DrawableTrimesh<> & m, int N, double l=0.05);
void non_uniform_tri(DrawableTrimesh<> & m,int N);
void make_domain(DrawableTrimesh<> & m, const double max_area, const double min_angle=20);
void make_triangulation(DrawableTrimesh<> &m, int mode, int N, double y_anisotropy=5);
void make_grid(DrawableTrimesh<> &m,int N);

#endif // TRIANGULATIONS_H
