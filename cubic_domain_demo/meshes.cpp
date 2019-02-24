#include <cinolib/meshes/meshes.h>
#include <cinolib/tetgen_wrap.h>
#include <cinolib/triangle_wrap.h>
#include "polygonsoup.h"

using namespace cinolib;

const int REDUNDANCY = 100;     // redundancy for Poisson disc sampling

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
// NON FUNZIONA - SERVE PROCEDURA DI SUDDIVISIONE MIDPOINT DI TETMESH - DA CHIEDERE A CINO

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
    /*double target_edge_length = 100.0/vert_per_edge;
    double max_tet_vol = 12*std::pow(target_edge_length,3.0)/6.0;       // ad hoc value to get ~vert_per_edge^3 tets
    std::string flags  = "q"+std::to_string(1.4)+"a" + std::to_string(max_tet_vol*z_anisotropy);
     //std::string flags  = "S" + std::to_string(vert_per_edge);
    tetgen_wrap({vec3d(0,0,0), vec3d(100,0,0), vec3d(100,0,100*z_anisotropy), vec3d(0,0,100*z_anisotropy),
                 vec3d(0,100,0), vec3d(100,100,0), vec3d(100,100,100*z_anisotropy), vec3d(0,100,100*z_anisotropy)},
                {{0,1,2,3},{4,5,6,7},{0,1,5,4},{3,2,6,7},{1,2,6,5},{0,3,7,4}},{},flags.c_str(),m);
    for (int i=0;i<m.num_verts();i++) m.vert(i).z() /= z_anisotropy;
    m.update_bbox();
    m.updateGL();*/
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

double my_rand()
{
    const double epsilon = 0.01;
    double r;
    do r = (double)rand()/RAND_MAX; while(r<epsilon || r>1.0-epsilon);
    return r;
}

void non_unif_Poisson_refine(int n, std::vector<vec3d> & buf, int frozen)
// refine samples contained in buf to return <= n samples
// with non-uniform density inv. prop. to 3D distance of point from origin
{
    double scale = 33;
    double r, rj, d;
    std::vector<bool> flagbuf(buf.size(),true);

    for (uint i=0;i<buf.size();i++)
    {
        if (!flagbuf[i]) continue;
        r = buf[i].length_squared()/scale;
        for (uint j=frozen;j<buf.size();j++)
        {
            if (i==j) continue;
            if (!flagbuf[j]) continue;
            rj = buf[j].length_squared()/scale;
            vec3d diff = buf[i]-buf[j];
            d = diff.length_squared();
            if (d<std::min(r,rj)) flagbuf[j]=false;
        }
    }

    int count = 0;
    for (uint i=0;i<flagbuf.size();i++)
    {
        if (flagbuf[i]) {buf[count]=buf[i]; count++;}
        if (count>=n) break;
    }
    buf.resize(count);
    std::cout << "Refined " << count << " points out of desired " << n << std::endl;
}

void Poisson1D(int N, vec3d v0, vec3d v1, int s, std::vector<vec3d> & e)
// s non-unif Poisson sample points along segment v0 v1; return points in e
{
    e = {v0, v1};
    if (v0.x()!=v1.x()) for (int i=0;i<s;i++) {e.push_back(vec3d(my_rand(),v0.y(),v0.z()));}
    if (v0.y()!=v1.y()) for (int i=0;i<s;i++) {e.push_back(vec3d(v0.x(),my_rand(),v0.z()));}
    if (v0.z()!=v1.z()) for (int i=0;i<s;i++) {e.push_back(vec3d(v0.x(),v0.y(),my_rand()));}
    non_unif_Poisson_refine(N, e, 2);
}

void make_hull(uint N, std::vector<vec3d> & points, std::vector<std::vector<uint>> & faces)
// sample points on the surface of the unit cube and make a mesh of the surface
{
    points={vec3d(0,0,0), vec3d(1,0,0), vec3d(1,0,1), vec3d(0,0,1),
            vec3d(0,1,0), vec3d(1,1,0), vec3d(1,1,1), vec3d(0,1,1)};

    std::vector<vec3d>  ex00, ex10, ex01, ex11, e0y0, e1y0, e0y1, e1y1, e00z, e10z, e01z, e11z;
    std::vector<vec3d>  fxy0,fxy1, fx0z, fx1z, f0yz, f1yz;

    int edgesamples = REDUNDANCY*N;
    int facesamples = REDUNDANCY*N*N;

    std::string flags;

    // sample 12 edges of cube
    Poisson1D(N,vec3d(0,0,0), vec3d(1,0,0),edgesamples,ex00);
    Poisson1D(N,vec3d(0,1,0), vec3d(1,1,0),edgesamples,ex10);
    Poisson1D(N,vec3d(0,0,1), vec3d(1,0,1),edgesamples,ex01);
    Poisson1D(N,vec3d(0,1,1), vec3d(1,1,1),edgesamples,ex11);
    Poisson1D(N,vec3d(0,0,0), vec3d(0,1,0),edgesamples,e0y0);
    Poisson1D(N,vec3d(1,0,0), vec3d(1,1,0),edgesamples,e1y0);
    Poisson1D(N,vec3d(0,0,1), vec3d(0,1,1),edgesamples,e0y1);
    Poisson1D(N,vec3d(1,0,1), vec3d(1,1,1),edgesamples,e1y1);
    Poisson1D(N,vec3d(0,0,0), vec3d(0,0,1),edgesamples,e00z);
    Poisson1D(N,vec3d(1,0,0), vec3d(1,0,1),edgesamples,e10z);
    Poisson1D(N,vec3d(0,1,0), vec3d(0,1,1),edgesamples,e01z);
    Poisson1D(N,vec3d(1,1,0), vec3d(1,1,1),edgesamples,e11z);

    // sample & mesh 6 faces of cube
    int frozen;
    std::vector<uint> segs = {0,1,1,2,2,3,3,0};

std::cout << "fxy0 ======================" << std::endl;
    fxy0={vec3d(0,0,0), vec3d(0,1,0), vec3d(1,1,0), vec3d(1,0,0)};
    fxy0.insert(fxy0.end(),ex00.begin()+2,ex00.end());
    fxy0.insert(fxy0.end(),ex10.begin()+2,ex10.end());
    fxy0.insert(fxy0.end(),e0y0.begin()+2,e0y0.end());
    fxy0.insert(fxy0.end(),e1y0.begin()+2,e1y0.end());
    frozen = fxy0.size();
    for (int i=0;i<facesamples;i++) fxy0.push_back(vec3d(my_rand(),my_rand(),0));
    non_unif_Poisson_refine(N*N, fxy0, frozen);
    std::vector<vec3d> mvxy0;
    std::vector<std::vector<uint>> mfxy0;
    //triangle_wrap(fxy0,segs,{},flags,mvxy0,mfxy0);
    fxy0.clear();
    for (uint i=0;i<mvxy0.size();i++) fxy0.push_back(vec3d(mvxy0[i].x(),mvxy0[i].y(),0));
    for (uint i=0;i<mfxy0.size();i++) {uint h=mfxy0[i][0]; mfxy0[i][0]=mfxy0[i][1]; mfxy0[i][1]=h;}


std::cout << "fxy1 ======================" << std::endl;
    fxy1={vec3d(0,0,1), vec3d(1,0,1), vec3d(1,1,1), vec3d(0,1,1)};
    fxy1.insert(fxy1.end(),ex01.begin()+2,ex01.end());
    fxy1.insert(fxy1.end(),ex11.begin()+2,ex11.end());
    fxy1.insert(fxy1.end(),e0y1.begin()+2,e0y1.end());
    fxy1.insert(fxy1.end(),e1y1.begin()+2,e1y1.end());
    frozen = fxy1.size();
    for (int i=0;i<facesamples;i++) fxy1.push_back(vec3d(my_rand(),my_rand(),1));
    non_unif_Poisson_refine(N*N, fxy1, frozen);
    std::vector<vec3d> mvxy1;
    std::vector<std::vector<uint>> mfxy1;
    //triangle_wrap(fxy1,segs,{},flags,mvxy1,mfxy1);
    fxy1.clear();
    for (uint i=0;i<mvxy1.size();i++) fxy1.push_back(vec3d(mvxy1[i].x(),mvxy1[i].y(),1));

std::cout << "f0yz ======================" << std::endl;
    f0yz={vec3d(0,0,0), vec3d(0,0,1), vec3d(0,1,1), vec3d(0,1,0)};
    f0yz.insert(f0yz.end(),e0y0.begin()+2,e0y0.end());
    f0yz.insert(f0yz.end(),e0y1.begin()+2,e0y1.end());
    f0yz.insert(f0yz.end(),e00z.begin()+2,e00z.end());
    f0yz.insert(f0yz.end(),e01z.begin()+2,e01z.end());
    frozen = f0yz.size();
    for (int i=0;i<facesamples;i++) f0yz.push_back(vec3d(0,my_rand(),my_rand()));
    non_unif_Poisson_refine(N*N, f0yz, frozen);
    std::vector<vec3d> mv0yz;
    std::vector<std::vector<uint>> mf0yz;
    for (uint i=0;i<f0yz.size();i++) f0yz[i].x()=f0yz[i].z();
    //triangle_wrap(f0yz,segs,{},flags,mv0yz,mf0yz);
    f0yz.clear();
    for (uint i=0;i<mv0yz.size();i++) f0yz.push_back(vec3d(0,mv0yz[i].y(),mv0yz[i].x()));

std::cout << "f1yz ======================" << std::endl;
    f1yz={vec3d(1,0,0), vec3d(1,1,0), vec3d(1,1,1), vec3d(1,0,1)};
    f1yz.insert(f1yz.end(),e1y0.begin()+2,e1y0.end());
    f1yz.insert(f1yz.end(),e1y1.begin()+2,e1y1.end());
    f1yz.insert(f1yz.end(),e10z.begin()+2,e10z.end());
    f1yz.insert(f1yz.end(),e11z.begin()+2,e11z.end());
    frozen = f1yz.size();
    for (int i=0;i<facesamples;i++) f1yz.push_back(vec3d(1,my_rand(),my_rand()));
    non_unif_Poisson_refine(N*N, f1yz, frozen);
    std::vector<vec3d> mv1yz;
    std::vector<std::vector<uint>> mf1yz;
    for (uint i=0;i<f1yz.size();i++) f1yz[i].x()=f1yz[i].z();
    //triangle_wrap(f1yz,segs,{},flags,mv1yz,mf1yz);
    f1yz.clear();
    for (uint i=0;i<mv1yz.size();i++) f1yz.push_back(vec3d(1,mv1yz[i].y(),mv1yz[i].x()));
    for (uint i=0;i<mf1yz.size();i++) {uint h=mf1yz[i][0]; mf1yz[i][0]=mf1yz[i][1]; mf1yz[i][1]=h;}

std::cout << "fx0z ======================" << std::endl;
    fx0z={vec3d(0,0,0), vec3d(1,0,0), vec3d(1,0,1), vec3d(0,0,1)};
    fx0z.insert(fx0z.end(),ex00.begin()+2,ex00.end());
    fx0z.insert(fx0z.end(),ex01.begin()+2,ex01.end());
    fx0z.insert(fx0z.end(),e00z.begin()+2,e00z.end());
    fx0z.insert(fx0z.end(),e10z.begin()+2,e10z.end());
    frozen = fx0z.size();
    for (int i=0;i<facesamples;i++) fx0z.push_back(vec3d(my_rand(),0,my_rand()));
    non_unif_Poisson_refine(N*N, fx0z, frozen);
    std::vector<vec3d> mvx0z;
    std::vector<std::vector<uint>> mfx0z;
    for (uint i=0;i<fx0z.size();i++) fx0z[i].y()=fx0z[i].z();
    //triangle_wrap(fx0z,segs,{},flags,mvx0z,mfx0z);
    fx0z.clear();
    for (uint i=0;i<mvx0z.size();i++) fx0z.push_back(vec3d(mvx0z[i].x(),0,mvx0z[i].y()));

std::cout << "fx1z ======================" << std::endl;
    fx1z={vec3d(0,1,0), vec3d(0,1,1), vec3d(1,1,1), vec3d(1,1,0)};
    fx1z.insert(fx1z.end(),ex11.begin()+2,ex11.end());
    fx1z.insert(fx1z.end(),ex10.begin()+2,ex10.end());
    fx1z.insert(fx1z.end(),e01z.begin()+2,e01z.end());
    fx1z.insert(fx1z.end(),e11z.begin()+2,e11z.end());
    frozen = fx1z.size();
    for (int i=0;i<facesamples;i++) fx1z.push_back(vec3d(my_rand(),1,my_rand()));
    non_unif_Poisson_refine(N*N, fx1z, frozen);
    std::vector<vec3d> mvx1z;
    std::vector<std::vector<uint>> mfx1z;
    for (uint i=0;i<fx1z.size();i++) fx1z[i].y()=fx1z[i].z();
    //triangle_wrap(fx1z,segs,{},flags,mvx1z,mfx1z);
    fx1z.clear();
    for (uint i=0;i<mvx1z.size();i++) fx1z.push_back(vec3d(mvx1z[i].x(),1,mvx1z[i].y()));

    // assemble hull
    polygonsoup s_hull;
    std::vector<vec3d> fbuf;

    for (uint i=0;i<mfxy0.size();i++) {
        fbuf.clear();
        for (uint j=0;j<mfxy0[i].size();j++) fbuf.push_back(fxy0[mfxy0[i][j]]);
        s_hull.add_polygon(fbuf);
    }
    for (uint i=0;i<mfxy1.size();i++) {
        fbuf.clear();
        for (uint j=0;j<mfxy1[i].size();j++) fbuf.push_back(fxy1[mfxy1[i][j]]);
        s_hull.add_polygon(fbuf);
    }
    for (uint i=0;i<mf0yz.size();i++) {
        fbuf.clear();
        for (uint j=0;j<mf0yz[i].size();j++) fbuf.push_back(f0yz[mf0yz[i][j]]);
        s_hull.add_polygon(fbuf);
    }
    for (uint i=0;i<mf1yz.size();i++) {
        fbuf.clear();
        for (uint j=0;j<mf1yz[i].size();j++) fbuf.push_back(f1yz[mf1yz[i][j]]);
        s_hull.add_polygon(fbuf);
    }
    for (uint i=0;i<mfx0z.size();i++) {
        fbuf.clear();
        for (uint j=0;j<mfx0z[i].size();j++) fbuf.push_back(fx0z[mfx0z[i][j]]);
        s_hull.add_polygon(fbuf);
    }
    for (uint i=0;i<mfx1z.size();i++) {
        fbuf.clear();
        for (uint j=0;j<mfx1z[i].size();j++) fbuf.push_back(fx1z[mfx1z[i][j]]);
        s_hull.add_polygon(fbuf);
    }

    DrawableTrimesh<> hull;
    s_hull.make_trimesh(hull);
//    hull.save("hull.obj");

    // return points and faces
    points = hull.vector_verts();
    faces = hull.vector_polys();
}


void non_unif_Poisson3D(int n, std::vector<vec3d> & buf)
// non-uniform Poisson disc sampling of n points in (0,1)x(0,1)x(0,1)
// with non-uniform density inv. prop. to distance of point from origin
// sampled points are ADDED to buf (no reset)
{
    int nsamples = REDUNDANCY*n;
    int frozen=buf.size();
    for (int i=0;i<nsamples;i++) buf.push_back(vec3d(my_rand(),my_rand(),my_rand()));
    non_unif_Poisson_refine(n, buf, frozen);
}

void non_uniform_tet(DrawableTetmesh<> & m,int N)
{
    std::vector<vec3d> points;
    std::vector<std::vector<uint>> faces={};
    std::string flags;

    make_hull(N,points,faces);
    non_unif_Poisson3D(N*N*N,points);
    tetgen_wrap(points,faces,{},flags.c_str(),m);
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
        non_uniform_tet(m,N);
        break;
    case 3:
        anisotropic_tet(m,N,anisotropy);
        break;
    }
}
