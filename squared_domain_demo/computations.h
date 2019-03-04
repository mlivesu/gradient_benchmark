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

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
bool poly_has_vert_on_boundary (const DrawableTrimesh<> &m, uint pid);
double relative_error(const vec3d v1, const vec3d v2, const int mode);
double absolute_error(const vec3d v1, const vec3d v2, const int mode);
double scale_function(double x, double md, double Md);
bool vector_contains_value(std::vector<uint> &v, const uint value);
double find_max_norm (const DrawableVectorField &V);
void find_max_min_values (const ScalarField f, double &max, double &min);
ScalarField heat_map_normalization(const ScalarField &f, double min, double max, double sat_neg=0, double sat_pos=1000);
DrawableVectorField compute_field(DrawableTrimesh<> &m, ScalarField & f, const int mode=0,bool handle_boundary=false);
ScalarField get_scalar_field(const DrawableTrimesh<> &m,const double a, const double b,const double c,const int mode=0);
DrawableVectorField  compute_ground_truth (const DrawableTrimesh<> &m,const double a, const double b,const double c,const int mode,const int method);
ScalarField estimate_error(const DrawableVectorField & GT, const DrawableVectorField & V, const DrawableTrimesh<> & m, int mode,int method,int relative,const int type_of_vertices);
std::vector<double> dual_error(const DrawableTrimesh<> &m, const DrawableVectorField &GT , const DrawableVectorField &V, const int mode, const int relative);
vec3d barycentric_coordinates(const vec3d &A,const vec3d &B, const vec3d &C, const vec3d &P);
void bring_the_field_inside(const DrawableTrimesh<> &m, DrawableTrimesh<> &m_grid, DrawableVectorField &V, DrawableVectorField &W, const int method);
Eigen::VectorXd fit_with_quadrics(const DrawableTrimesh<> & m, uint i, ScalarField &f);

//VERTEX-BASED METHODS
Eigen::SparseMatrix<double> build_matrix_for_AGS(const DrawableTrimesh<> & m);
void build_matrix_for_LSDD(DrawableTrimesh<> &m, std::vector<Eigen::ColPivHouseholderQR<Eigen::MatrixXd> > &MFact, std::vector<Eigen::MatrixXd> &RHS, std::vector<std::vector<uint> > &nbrs);
void build_matrix_for_LR(DrawableTrimesh<> &m, std::vector<Eigen::ColPivHouseholderQR<Eigen::MatrixXd> > &M, std::vector<Eigen::MatrixXd> &RHS, std::vector<std::vector<uint> > &nbrs);
void solve_for_LSDD(DrawableTrimesh<> &m, DrawableVectorField &V, std::vector<Eigen::ColPivHouseholderQR<Eigen::MatrixXd> > &M, std::vector<Eigen::MatrixXd> &RHS,  ScalarField & f, std::vector<std::vector<uint> > &nbrs);
void solve_for_LR(DrawableTrimesh<> &m, DrawableVectorField &V, std::vector<Eigen::ColPivHouseholderQR<Eigen::MatrixXd> > &M, std::vector<Eigen::MatrixXd> RHS,  ScalarField & f, std::vector<std::vector<uint> > &nbrs);
#endif // COMPUTATIONS_H
