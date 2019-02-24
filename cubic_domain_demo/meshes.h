#ifndef MESHES_H
#define MESHES_H

#include <cinolib/meshes/meshes.h>

using namespace cinolib;

void make_tet_mesh(DrawableTetmesh<> & m, int vert_per_edge, int mode,double anisotropy=3);

#endif // MESHES_H
