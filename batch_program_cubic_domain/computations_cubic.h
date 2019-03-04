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

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
bool vector_contains_value(const std::vector<uint> &v,const uint value);
bool poly_has_vert_on_srf (const DrawableTetmesh<> &m, uint pid);
double relative_error(const vec3d v1, const vec3d v2, const int mode);
double absolute_error(const vec3d v1, const vec3d v2, const int mode);
double scale_function(double x, double md, double Md,const double scale_factor);


uint closest_vertex (const vec3d & p, const DrawableTrimesh<> *m);
DrawableVectorField arrows_normalization(const DrawableTrimesh<> &m, const DrawableVectorField &V, const int mode, const int scale_factor);
double find_max_norm(const DrawableVectorField &V);
void find_max_min_values (const ScalarField f, double &max, double &min);
ScalarField heat_map_normalization(const ScalarField &f,double min,double max,double scale_factor);


DrawableVectorField compute_field(const DrawableTetmesh<> &m, ScalarField & f, const int method,int weight=0);
ScalarField  get_scalar_field(const DrawableTetmesh<> &m, const double a, const double b, const double c, const int method, int noise, double k);

DrawableVectorField compute_quadratic_regression(const DrawableTetmesh<> &m, const ScalarField & f);
DrawableVectorField compute_quadratic_regression_with_centroids(const DrawableTetmesh<> &m, const ScalarField & f);

DrawableVectorField compute_FEM(const DrawableTetmesh<> &m, const ScalarField & f);
DrawableVectorField compute_FEM_with_centroids(const DrawableTetmesh<> &m, const ScalarField & f);
//Eigen::SparseMatrix<double> GG_gradient_matrix(const DrawableTetmesh<> &m, int method);
DrawableVectorField from_p2v(DrawableVectorField &W, const DrawableTetmesh<> &m,int weight=0);

double solid_angle(const DrawableTetmesh<> &m,uint &pid, uint &vid);

DrawableVectorField compute_ground_truth(const DrawableTetmesh<> &m, const int method, const double a, const double b, const double c=10, const int mode=2);

double estimate_MSE(const DrawableVectorField &GT, const DrawableVectorField &V, const DrawableTetmesh<> &grid, const int method, const int type_of_vertices, int mode);
double dual_error(const DrawableTetmesh<> &m, const DrawableVectorField &GT , const DrawableVectorField &V, const int mode);
std::vector<double> circum_radius(const DrawableTrimesh<> &m);
std::vector<double> error_for_hist(const DrawableVectorField &GT, const DrawableVectorField &V, const DrawableTetmesh<> &m, const int method, const int type_of_vertices, int mode);
std::vector<double> in_radius(const DrawableTrimesh<> &m);
std::vector<double> radii_ratio(const DrawableTrimesh<> &m, const std::vector<double> &R, std::vector<double> &r );
double correlation_coefficient(const std::vector<double> &X, const std::vector<double> &Y);
double average_neighborhood_area(const DrawableTrimesh<> &m);
void bring_the_field_inside(const DrawableTetmesh<> &m, DrawableTetmesh<> &m_grid, DrawableVectorField &V, DrawableVectorField &W, const int method);
std::vector<double> barycentric_coordinates(const vec3d &A,const vec3d &B, const vec3d &C, const vec3d &D, const vec3d &P);
DrawableVectorField compute_PCE(const DrawableTetmesh<> & m, const ScalarField &f);
DrawableVectorField compute_AGS(const DrawableTetmesh<> & m, ScalarField &f, int weight=0);
std::vector<double> compute_for_hist(DrawableTetmesh<> &m, DrawableTetmesh<> &grid,const int method,const double a, const double b,const int tri_type, const int N, ScalarField &F,const double anisotropy,int without_boundary,int mode,int weight);
Eigen::SparseMatrix<double> build_matrix_for_AGS(const DrawableTetmesh<> & m);
void build_matrix_for_LSDD(DrawableTetmesh<> &m, std::vector<Eigen::ColPivHouseholderQR<Eigen::MatrixXd> > &MFact, std::vector<Eigen::MatrixXd> &RHS, std::vector<std::vector<uint> > &nbrs);
void build_matrix_for_LR(DrawableTetmesh<> &m, std::vector<Eigen::ColPivHouseholderQR<Eigen::MatrixXd> > &M, std::vector<Eigen::MatrixXd> &RHS,std::vector<std::vector<uint> > &nbrs);
void solve_for_LSDD(DrawableTetmesh<> &m, DrawableVectorField &V, std::vector<Eigen::ColPivHouseholderQR<Eigen::MatrixXd> > &M, std::vector<Eigen::MatrixXd> &RHS, ScalarField & f, std::vector<std::vector<uint> > &nbrs);
void solve_for_LR(DrawableTetmesh<> &m,DrawableVectorField &V, std::vector<Eigen::ColPivHouseholderQR<Eigen::MatrixXd> > &M, std::vector<Eigen::MatrixXd> &RHS, ScalarField &f, std::vector<std::vector<uint> > &nbrs);
#endif // COMPUTATIONS_H
