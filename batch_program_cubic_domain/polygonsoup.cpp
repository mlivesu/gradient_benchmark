#include "polygonsoup.h"

class vertex_sorter
{
public:
    vec3d coords;
    uint face;
    uint in_face;

    vertex_sorter(vec3d v,uint f, uint i) { coords = v; face = f; in_face = i; }
};
bool compare(vertex_sorter a, vertex_sorter b) { return(a.coords < b.coords); }

polygonsoup::polygonsoup() { polys.clear(); }

polygonsoup::polygonsoup(DrawableTrimesh<> &m)
{
    polys.clear();
    for (uint i=0;i<m.num_polys();i++)
    {
        polys.push_back(m.poly_vlist(i));
    }
}

void polygonsoup::add_polygon(std::vector<vec3d> p)
{
    polys.push_back(p);
}

bool same_point(vec3d v, vec3d w)
{
    return (v.x()==w.x() && v.y()==w.y() && v.z()==w.z());
}

void polygonsoup::make_trimesh(DrawableTrimesh<> &m)
{
   std::vector<vertex_sorter> buf;
   std::vector<vec3d> mverts;
   std::vector<std::vector<uint>> mfaces;
   m.clear();

   // allocate faces for indexed mesh
   for (uint i=0;i<polys.size();i++)  mfaces.push_back(std::vector<uint>(polys[i].size()));

   // sort vertices to get contiguous duplicates
   for(uint i=0;i<polys.size();i++)
       for(uint j=0;j<polys[i].size();j++)
           buf.push_back(vertex_sorter(polys[i][j],i,j));
   std::sort(buf.begin(),buf.end(),compare);

   // fill-in vertices and faces of indexed mesh
   int curi = 0;                                    // current vertex index
   uint i=0;
   while (i<buf.size())
   {
       mverts.push_back(buf[i].coords);             // add new vertex
       mfaces[buf[i].face][buf[i].in_face]=curi;    // set incidence in face
       i++;
       while(i<buf.size() && same_point(buf[i].coords,mverts[curi])) // set incidence for duplicates
       {
           mfaces[buf[i].face][buf[i].in_face]=curi;
           i++;
       }
       curi++;
   }

   // build mesh
   m = DrawableTrimesh<>(mverts,mfaces);
}
