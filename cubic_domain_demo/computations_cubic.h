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

    DrawableTetmesh<>                m;
    DrawableTetmesh<>                m_grid;
    DrawablePolyhedralmesh<>         dual_m;
    DrawablePolyhedralmesh<>         dual_m_grid;

    DrawableVectorField              V;
    DrawableVectorField              V_norm;
    DrawableVectorField              V_grid;

    DrawableVectorField               GT;
    DrawableVectorField               GT_grid;
    DrawableVectorField               GT_norm;

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
bool vector_contains_value(std::vector<uint> &v, const uint value);
bool poly_has_vert_on_srf (const DrawableTetmesh<> &m, uint pid);
double relative_error(const vec3d v1, const vec3d v2, const int mode);
double absolute_error(const vec3d v1, const vec3d v2, const int mode);
double scale_function(double x, double md, double Md);
double range(const std::vector<double> v);

uint closest_vertex (const vec3d & p, const DrawablePolyhedralmesh<> *m);
uint closest_vertex (const vec3d & p, const DrawableTetmesh<> *m);
DrawableVectorField arrows_normalization(const DrawableTrimesh<> &m, const DrawableVectorField &V, const int mode, const int scale_factor);
double find_max_norm(const DrawableVectorField &V);
void find_max_min_values (const ScalarField f, double &max, double &min);
ScalarField heat_map_normalization(const ScalarField &f, double min, double max, double sat_neg=0, double sat_pos=1000);


DrawableVectorField compute_field(DrawableTetmesh<> &m, ScalarField & f, const int method);
ScalarField get_scalar_field(const DrawableTetmesh<> &m, const double a, const double b, const double c, const int method, int noise, double k);




DrawableVectorField compute_ground_truth(const DrawableTetmesh<> &m, const double a, const double b, const double c, const int mode, const int method);
ScalarField estimate_error(const DrawableVectorField & GT, const DrawableVectorField & V, const DrawableTetmesh<> &m, const int mode, const int method, const int relative, const int type_of_vertices);

double estimate_MSE(const DrawableVectorField &GT, const DrawableVectorField &V, const DrawableTetmesh<> &m, const int method, const int type_of_vertices, int mode);
std::vector<double> dual_error(const DrawableTetmesh<> &m, const DrawableVectorField &GT , const DrawableVectorField &V, const int mode, const int relative);
void bring_the_field_inside(const DrawableTetmesh<> &m, DrawableTetmesh<> &m_grid, DrawableVectorField &V, DrawableVectorField &W, const int method);
std::vector<double> barycentric_coordinates(const vec3d &A,const vec3d &B, const vec3d &C, const vec3d &D, const vec3d &P);

//VERTEX-BASED METHODS
Eigen::SparseMatrix<double> build_matrix_for_AGS(const DrawableTetmesh<> & m);
void build_matrix_for_LSDD(DrawableTetmesh<> &m, std::vector<Eigen::ColPivHouseholderQR<Eigen::MatrixXd> > &MFact, std::vector<Eigen::MatrixXd> &RHS, std::vector<std::vector<uint> > &nbrs);
void build_matrix_for_LR(DrawableTetmesh<> &m, std::vector<Eigen::ColPivHouseholderQR<Eigen::MatrixXd> > &M, std::vector<Eigen::MatrixXd> &RHS,std::vector<std::vector<uint> > &nbrs);
void solve_for_LSDD(DrawableTetmesh<> &m, DrawableVectorField &V, std::vector<Eigen::ColPivHouseholderQR<Eigen::MatrixXd> > &M, std::vector<Eigen::MatrixXd> &RHS, ScalarField & f, std::vector<std::vector<uint> > &nbrs);
void solve_for_LR(DrawableTetmesh<> &m,DrawableVectorField &V, std::vector<Eigen::ColPivHouseholderQR<Eigen::MatrixXd> > &M, std::vector<Eigen::MatrixXd> &RHS, ScalarField &f, std::vector<std::vector<uint> > &nbrs);

#endif // COMPUTATIONS_H
