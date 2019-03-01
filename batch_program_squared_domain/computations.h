#ifndef COMPUTATIONS_H
#define COMPUTATIONS_H


#include <cinolib/meshes/meshes.h>
#include <cinolib/scalar_field.h>
#include <cinolib/drawable_vector_field.h>
#include <cinolib/scalar_field.h>
#include <Eigen/Sparse>
#include <cinolib/geometry/vec3.h>
#include <cinolib/dual_mesh.h>
#include <cinolib/profiler.h>
#include <cinolib/cino_inline.h>
#include <cinolib/gradient.h>


using namespace cinolib;
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

typedef struct{

    DrawableTrimesh<>                m;
    DrawablePolygonmesh<>            dual_m;
    DrawableVectorField   V;
    DrawableVectorField   V_norm;
    DrawableVectorField   GT;
    DrawableVectorField   GT_norm;

    std::vector<int>                 rank;


    ScalarField                      f;

    ScalarField                      f_norm;

    ScalarField                      err;


    ScalarField                      err_norm;



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
double scale_function(double x, const double scale_factor);
uint closest_vertex (const vec3d & p, const DrawableTrimesh<> *m);
DrawableVectorField arrows_normalization(const DrawableTrimesh<> &m, const DrawableVectorField &V, const int mode, const int scale_factor);
void find_max_min_values (const DrawableVectorField &V, vec3d &max, vec3d &min);
void find_max_min_values (const ScalarField f, double &max, double &min);
ScalarField heat_map_normalization(const ScalarField &f, double scale_factor);
DrawableVectorField compute_field(DrawableTrimesh<> &m, ScalarField & f, const int method,int weight=0);
std::vector<double> set_coefficients(const std::vector<vec3d> coords, const std::vector<double> weights,const std::vector<double> Aij);
double sum_up_value(const std::vector<double> coords, const std::vector<double> weights,const std::vector<double> Aij);
DrawableVectorField compute_quadratic_regression(const DrawableTrimesh<> &m, const ScalarField & f);
DrawableVectorField compute_quadratic_regression_centroids(DrawableTrimesh<> &m, const ScalarField & f);
DrawableVectorField compute_FEM(const DrawableTrimesh<> &m, const ScalarField & f);
ScalarField get_scalar_field(const DrawableTrimesh<> &m, const double a, const double b, const int method, int noise, double k);
DrawableVectorField  compute_ground_truth (const DrawableTrimesh<> &m,const int method,const double a, const double b,const double c=10,const int mode=2);
DrawableVectorField from_f2v(DrawableVectorField &W, DrawableTrimesh<> &m,int weight=0);
double estimate_MSE(const DrawableVectorField & GT, const DrawableVectorField & V, const DrawableTrimesh<> & m, const int method, const int type_of_vertices, int mode);
double dual_error(const DrawableTrimesh<> &m, const DrawableVectorField &GT , const DrawableVectorField &V, const int mode);
std::vector<double> circum_radius(const DrawableTrimesh<> &m);
std::vector<double> in_radius(const DrawableTrimesh<> &m);
std::vector<double> radii_ratio(const DrawableTrimesh<> &m, const std::vector<double> &R, std::vector<double> &r );
double correlation_coefficient(const std::vector<double> &X, const std::vector<double> &Y);
double average_neighborhood_area(const DrawableTrimesh<> &m);
std::vector<double> error_for_hist(const DrawableVectorField &GT, const DrawableVectorField &V, const DrawableTrimesh<> &m, const int method, const int type_of_vertices, int mode);
void bring_the_field_inside(const DrawableTrimesh<> &m, DrawableTrimesh<> &m_grid, DrawableVectorField &V, DrawableVectorField &W, const int method);
vec3d barycentric_coordinates(const vec3d &A,const vec3d &B, const vec3d &C, const vec3d &P);
DrawableVectorField compute_PCE(const DrawableTrimesh<> & m, ScalarField &f);
DrawableVectorField compute_AGS(const DrawableTrimesh<> & m, ScalarField &f, int weight=0);
DrawableVectorField GG_on_verts(const DrawableTrimesh<> & m, ScalarField &f);
std::vector<double> compute_for_hist(DrawableTrimesh<> &m, DrawableTrimesh<> &grid,const int method,const double a, const double b,const int tri_type, const int N, ScalarField &F,const double anisotropy,int without_boundary,int mode,int weight);
void build_matrix_for_LSDD(DrawableTrimesh<> &m, std::vector<Eigen::ColPivHouseholderQR<Eigen::Matrix2d> > &MFact, std::vector<Eigen::MatrixXd> &RHS, std::vector<std::vector<uint> > &nbrs);
void build_matrix_for_LR(DrawableTrimesh<> &m, std::vector<Eigen::ColPivHouseholderQR<Eigen::MatrixXd> > &M, std::vector<Eigen::MatrixXd> &RHS,std::vector<std::vector<uint> > &nbrs);
void factorize(std::vector<Eigen::MatrixXd> M);
void solve_for_LSDD(DrawableTrimesh<> &m, DrawableVectorField &V, std::vector<Eigen::ColPivHouseholderQR<Eigen::Matrix2d> > &M, std::vector<Eigen::MatrixXd> &RHS, ScalarField & f, std::vector<std::vector<uint> > &nbrs, std::chrono::duration<double> time_precom, std::chrono::duration<double> time_estimation);
void solve_for_LR(DrawableTrimesh<> &m,DrawableVectorField &V, std::vector<Eigen::ColPivHouseholderQR<Eigen::MatrixXd> > &M, std::vector<Eigen::MatrixXd> RHS, ScalarField &f, std::vector<std::vector<uint> > &nbrs, std::chrono::duration<double> time_precom, std::chrono::duration<double> time_estimation);
#endif // COMPUTATIONS_H
