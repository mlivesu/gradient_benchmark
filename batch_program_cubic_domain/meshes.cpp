#include <cinolib/meshes/meshes.h>
#include <cinolib/tetgen_wrap.h>
#include <cinolib/triangle_wrap.h>
#include "polygonsoup.h"

using namespace cinolib;

const int REDUNDANCY = 100;     // redundancy for Poisson disc sampling

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


void regular_tet(DrawableTetmesh<> & m, const int vert_per_edge)
{
    std::vector<vec3d> points;
    std::vector<std::vector<uint>> dummy;

    int n=vert_per_edge+1;
    double step = 100.0/n;

    for (uint i=0;i<=n;i++)
        for (uint j=0;j<=n;j++)
            for (uint k=0;k<=n;k++)
                points.push_back(vec3d(i*step,j*step,k*step));

    /*std::string flags  = "pc";
    tetgen_wrap(points,dummy,{},flags.c_str(),m);
    m.updateGL();*/
    double target_edge_length = 100.0/vert_per_edge;
    double radius_edge_ratio=sqrt(3/8);
    double min_angle=atan(2*sqrt(2));
    double max_tet_vol=std::pow(target_edge_length,3.0)/(6*sqrt(2));
    //std::string flags  = "q"+std::to_string(radius_edge_ratio)+"/"+ std::to_string(min_angle);//+"a"+std::to_string(max_tet_vol);
    std::string flags  = "q"+std::to_string(radius_edge_ratio)+"/"+ std::to_string(min_angle)+"c";
    tetgen_wrap(points,dummy,{},flags.c_str(),m);
    m.updateGL();

}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void uniform_tet(DrawableTetmesh<> & m, const int vert_per_edge)
{
    double target_edge_length = 100.0/vert_per_edge;
    double max_tet_vol = 12*std::pow(target_edge_length,3.0)/6.0;       // ad hoc value to get ~vert_per_edge^3 tets
    std::string flags  = "q1.4a" + std::to_string(max_tet_vol);
    tetgen_wrap({vec3d(0,0,0), vec3d(100,0,0), vec3d(100,0,100), vec3d(0,0,100),
                 vec3d(0,100,0), vec3d(100,100,0), vec3d(100,100,100), vec3d(0,100,100)},
                {{0,1,2,3},{4,5,6,7},{0,1,5,4},{3,2,6,7},{1,2,6,5},{0,3,7,4}},{},flags.c_str(),m);
    m.updateGL();
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void anisotropic_tet(DrawableTetmesh<> & m, const int vert_per_edge, const double z_anisotropy)
{

    DrawableTetmesh<> aux;
    uniform_tet(aux,vert_per_edge);
    int Nv=aux.num_verts();

    double target_edge_length = 100.0/vert_per_edge;
    double max_tet_vol = 12*std::pow(target_edge_length,3.0)/6.0;       // ad hoc value to get ~vert_per_edge^3 tets

    std::string flags  = "S" + std::to_string(Nv-8)+"a" + std::to_string(max_tet_vol*z_anisotropy);
    tetgen_wrap({vec3d(0,0,0), vec3d(100,0,0), vec3d(100,0,100*z_anisotropy), vec3d(0,0,100*z_anisotropy),
                 vec3d(0,100,0), vec3d(100,100,0), vec3d(100,100,100*z_anisotropy), vec3d(0,100,100*z_anisotropy)},
                {{0,1,2,3},{4,5,6,7},{0,1,5,4},{3,2,6,7},{1,2,6,5},{0,3,7,4}},{},flags.c_str(),m);
    for (int i=0;i<m.num_verts();i++) m.vert(i).z() /= z_anisotropy;
    m.update_bbox();
    m.updateGL();


}




//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void make_tet_mesh(DrawableTetmesh<> & m, int N, int mode, double anisotropy)
{
    switch(mode)
    {
    case 0:
        regular_tet(m,N);
        break;
    case 1:
        uniform_tet(m,N);
        break;
    case 2:
        anisotropic_tet(m,N,anisotropy);
        break;

    }
}
