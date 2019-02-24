#ifndef POLYGONSOUP_H
#define POLYGONSOUP_H

#include <vector>
#include <cinolib/geometry/vec3.h>
#include <cinolib/meshes/meshes.h>

using namespace cinolib;

class polygonsoup
{
public:
    polygonsoup();

    polygonsoup(DrawableTrimesh<> & m);

    void add_polygon(std::vector<vec3d> p);

    void make_trimesh(DrawableTrimesh<> & m);

private:
    std::vector<std::vector<vec3d>> polys;
};

#endif // POLYGONSOUP_H
