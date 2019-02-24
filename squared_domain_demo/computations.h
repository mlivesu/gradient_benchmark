#ifndef COMPUTATIONS_H
#define COMPUTATIONS_H


#include <cinolib/meshes/meshes.h>
#include <cinolib/scalar_field.h>
#include <cinolib/drawable_vector_field.h>
#include <cinolib/scalar_field.h>
#include <Eigen/Sparse>
#include <cinolib/geometry/vec3.h>
#include <cinolib/dual_mesh.h>



using namespace cinolib;

enum contributes
{
    generic,
    right,
    left
};
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

typedef struct{

    DrawableTrimesh<>                m;
    DrawableTrimesh<>                m_grid;

    DrawablePolygonmesh<>            dual_m;

    DrawableVectorField              V;
    DrawableVectorField              V_norm;
    DrawableVectorField              V_grid;
    DrawableVectorField              GT;
    DrawableVectorField              GT_grid;
    DrawableVectorField              GT_norm;

    std::vector<int>                 rank;

    ScalarField                      f;

    ScalarField                      f_norm;

    ScalarField                      err;


    ScalarField                      err_norm;

    double                           max;

    double                           min;

    Eigen::SparseMatrix<double>      G;
}
COMPUTATIONS;

typedef struct{
    vec3d p;
    uint real;
    }
imaginary_vertex;
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
bool poly_has_vert_on_srf (const DrawableTrimesh<> &m, uint pid);
double relative_error(const vec3d v1, const vec3d v2, const int mode);
double absolute_error(const vec3d v1, const vec3d v2, const int mode);
double scale_function(double x, double md, double Md);
std::vector<imaginary_vertex> nbr_for_boundaries(const DrawableTrimesh<> &m, uint vid);
DrawableVectorField arrows_normalization(const DrawableTrimesh<> &m, const DrawableVectorField &V, const int mode, const int scale_factor);
double find_max_norm (const DrawableVectorField &V);
void find_max_min_values (const ScalarField f, double &max, double &min);
ScalarField heat_map_normalization(const ScalarField &f, double min, double max, double sat_neg=0, double sat_pos=1000);
DrawableVectorField compute_field(DrawableTrimesh<> &m, ScalarField & f, const int mode=0,bool handle_boundary=false);
std::vector<double> set_coefficients(const std::vector<vec3d> coords, const std::vector<double> weights,const std::vector<double> Aij);
double sum_up_value(const std::vector<double> coords, const std::vector<double> weights,const std::vector<double> Aij);
DrawableVectorField compute_quadratic_regression(DrawableTrimesh<> &m, const ScalarField & f,bool handle_boundaries=false);
DrawableVectorField compute_quadratic_regression_centroids(DrawableTrimesh<> &m, const ScalarField & f,bool handle_boundaries=false);
DrawableVectorField compute_FEM(DrawableTrimesh<> &m, const ScalarField & f, bool handle_bondary=false);
ScalarField get_scalar_field(const DrawableTrimesh<> &m,const double a, const double b,const double c,const int mode=0);
ScalarField scalar_field_with_boundaries(const DrawableTrimesh<> &m,const double a, const double b,const double c,const int mode,std::map<uint,int> &boundaries);
Eigen::SparseMatrix<double> GG_gradient_matrix(const DrawableTrimesh<> &m,int mode=0);
DrawableVectorField  compute_ground_truth (const DrawableTrimesh<> &m,const double a, const double b,const double c,const int mode,const int method);
DrawableVectorField from_f2v(DrawableVectorField &W, DrawableTrimesh<> &m);
ScalarField estimate_error(const DrawableVectorField & GT, const DrawableVectorField & V, const DrawableTrimesh<> & m, int mode,int method,int relative,const int type_of_vertices);
std::vector<double> dual_error(const DrawableTrimesh<> &m, const DrawableVectorField &GT , const DrawableVectorField &V, const int mode, const int relative);
vec3d barycentric_coordinates(const vec3d &A,const vec3d &B, const vec3d &C, const vec3d &P);
void bring_the_field_inside(const DrawableTrimesh<> &m, DrawableTrimesh<> &m_grid, DrawableVectorField &V, DrawableVectorField &W, const int method);
vec3d generic_contribute(const vec3d vi, const vec3d vj, const vec3d vh, const vec3d vk, const vec3d Njh, const vec3d Njk, int mode);
vec3d vi_contribute(const vec3d vi, const vec3d vj, const vec3d vh, const vec3d vk, const vec3d Njh, const vec3d Njk, int mode);
void pid_contribute(const std::vector<imaginary_vertex> &nbr, std::map<uint, int> boundaries, std::vector<std::pair<uint,vec3d>> &pesi, int &offset,vec3d &n);
DrawableVectorField compute_PCE(const DrawableTrimesh<> & m, ScalarField &f);
DrawableVectorField compute_AGS(const DrawableTrimesh<> & m, ScalarField &f);
void num_verts_on_boundary(const DrawableTrimesh<> & m, std::map<uint,int> &boundaries);
Eigen::VectorXd fit_with_quadrics(const DrawableTrimesh<> & m, uint i, ScalarField &f);
#endif // COMPUTATIONS_H
