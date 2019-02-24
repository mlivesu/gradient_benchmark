#ifndef PICKING_H
#define PICKING_H

#include <cinolib/gui/qt/glcanvas.h>
#include <cinolib/geometry/vec3.h>
#include <cinolib/meshes/meshes.h>
#include "computations.h"
using namespace cinolib;

void handle_mouse_press(cinolib::GLcanvas *c, QMouseEvent *e);
uint closest_vertex (const vec3d & p, const DrawablePolygonmesh<> *m);
uint closest_vertex (const vec3d & p, const DrawableTrimesh<> *m);
std::vector<uint> closest_vertices (const vec3d & p, const DrawableTrimesh<> *m);
#endif // PICKING_H
