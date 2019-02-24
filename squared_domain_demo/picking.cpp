#include "picking.h"
//comment to test githbub updating
using namespace cinolib;

extern COMPUTATIONS computations;
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
uint closest_vertex(const vec3d & p, const DrawablePolygonmesh<>  *m)
{
    std::vector<std::pair<double,uint>> closest;
    for(uint vid=0; vid<m->num_verts(); ++vid) closest.push_back(std::make_pair(m->vert(vid).dist(p),vid));
    std::sort(closest.begin(), closest.end());
    return closest.front().second;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
uint closest_vertex (const vec3d & p, const DrawableTrimesh<> *m)
{
    std::vector<std::pair<double,uint>> closest;
    for(uint vid=0; vid<m->num_verts(); ++vid) closest.push_back(std::make_pair(m->vert(vid).dist(p),vid));
    std::sort(closest.begin(), closest.end());
    return closest.front().second;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
std::vector<uint> closest_vertices(const vec3d & p, const DrawableTrimesh<> *m)
{
    std::vector<std::pair<double,uint>> closest;
    std::vector<uint> three_closest;
    for(uint vid=0; vid<m->num_verts(); ++vid) closest.push_back(std::make_pair(m->vert(vid).dist(p),vid));
    std::sort(closest.begin(), closest.end());
    for(std::vector<std::pair<double,uint>>::const_iterator iter=closest.begin();iter!=closest.end();++iter)
    {
        three_closest.push_back(iter->second);
    }


    return three_closest;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void handle_mouse_press(cinolib::GLcanvas *c, QMouseEvent *e)
{
    /*//WORK IN PROGRESS per implementarlo nelle 4 canvas

        if(e->modifiers() == Qt::ControlModifier)
        {

            vec3d p;
            vec2i click(e->x(), e->y());
            if (c->unproject(click, p)) // transform click in a 3d point
            {

            uint vid = closest_vertex(p,&computations.m);
            std::cout << computations.rank[vid] << std::endl;
        }
        }
        else if(e->modifiers() == Qt::ShiftModifier)
        {
            std::cout << "SHIFT+click: do something else!" << std::endl;
        }
        c->updateGL();*/

}
