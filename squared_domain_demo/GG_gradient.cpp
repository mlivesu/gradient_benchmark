/*********************************************************************************
*  Copyright(C) 2016: Marco Livesu                                               *
*  All rights reserved.                                                          *
*                                                                                *
*  This file is part of CinoLib                                                  *
*                                                                                *
*  CinoLib is dual-licensed:                                                     *
*                                                                                *
*   - For non-commercial use you can redistribute it and/or modify it under the  *
*     terms of the GNU General Public License as published by the Free Software  *
*     Foundation; either version 3 of the License, or (at your option) any later *
*     version.                                                                   *
*                                                                                *
*   - If you wish to use it as part of a commercial software, a proper agreement *
*     with the Author(s) must be reached, based on a proper licensing contract.  *
*                                                                                *
*  This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE       *
*  WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.     *
*                                                                                *
*  Author(s):                                                                    *
*                                                                                *
*     Marco Livesu (marco.livesu@gmail.com)                                      *
*     http://pers.ge.imati.cnr.it/livesu/                                        *
*                                                                                *
*     Italian National Research Council (CNR)                                    *
*     Institute for Applied Mathematics and Information Technologies (IMATI)     *
*     Via de Marini, 6                                                           *
*     16149 Genoa,                                                               *
*     Italy                                                                      *
**********************************************************************************/
#include "GG_gradient.h"
#include <cinolib/meshes/meshes.h>

typedef Eigen::Triplet<double> Entry;

using namespace cinolib;

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

template<class M, class V, class E, class P>
CINO_INLINE
Eigen::SparseMatrix<double> GG_gradient_matrix(const AbstractPolygonMesh<M,V,E,P> & m,int mode)
{
    if (mode==0)
    {
        Eigen::SparseMatrix<double> G(m.num_polys()*3, m.num_verts());
        std::vector<Entry> entries;
        
        for(uint pid=0; pid<m.num_polys(); ++pid)
        {
            double area = std::max(m.poly_area(pid), 1e-5) * 2.0; // (2 is the average term : two verts for each edge)
            vec3d n     = m.poly_data(pid).normal;
            
            for(uint off=0; off<m.verts_per_poly(pid); ++off)
            {
                uint  prev = m.poly_vert_id(pid,off);
                uint  curr = m.poly_vert_id(pid,(off+1)%m.verts_per_poly(pid));
                uint  next = m.poly_vert_id(pid,(off+2)%m.verts_per_poly(pid));
                vec3d u    = m.vert(next) - m.vert(curr);
                vec3d v    = m.vert(curr) - m.vert(prev);
                vec3d u_90 = u.cross(n); u_90.normalize();
                vec3d v_90 = v.cross(n); v_90.normalize();
                
                vec3d per_vert_sum_over_edge_normals = u_90 * u.length() + v_90 * v.length();
                per_vert_sum_over_edge_normals /= area;
                
                uint row = 3 * pid;
                entries.push_back(Entry(row, curr, per_vert_sum_over_edge_normals.x())); ++row;
                entries.push_back(Entry(row, curr, per_vert_sum_over_edge_normals.y())); ++row;
                entries.push_back(Entry(row, curr, per_vert_sum_over_edge_normals.z()));
            }
        }
        
        G.setFromTriplets(entries.begin(), entries.end());
        return G;

    }
    else if (mode==1){
        Eigen::SparseMatrix<double> G(m.num_verts()*3, m.num_verts());
        std::vector<Entry> entries;
        uint Nf=m.num_polys();
        uint Nv=m.num_verts();
        std::vector<uint> nbr;
        std::vector<uint> sx;
        std::vector<uint> dx;
        vec3d nij;
        vec3d nji;
        vec3d nik;
        vec3d nhi;
        vec3d Njh;
        vec3d Njk;
        vec3d value;
        vec3d vertex_contribute;
        vec3d first;
        vec3d second;
        vec3d third;

        double area=0.0;
        std::vector<vec3d> normals(Nf);
        std::vector<vec3d> vertex_normals(Nv);
        std::vector<vec3d> coords=m.vector_verts();
        int row=0;

        for (uint i=0; i<Nf;++i)
        {
            normals[i]=m.poly_data(i).normal;
            normals[i].normalize();
        }

        for(uint j=0;j<Nv;++j)
        {
            vertex_normals[j]=m.vert_data(j).normal;
            vertex_normals[j].normalize();
        }

        
        for (uint i=0;i<Nv;++i)
        {
            if(m.vert_is_boundary(i))
            {
                /*nbr=m.vert_ordered_vert_ring(i);
                area=0.0;
                for (int j=0;j<nbr.size();++j)
                {
                    int eid=m.edge_id(i,nbr[j]);
                    std::vector<uint>facce=m.adj_e2p(eid);
                    if(j==0)
                    {
                        vec3d n=normals[facce[0]];
                        nji=n.cross(((coords[i].operator +(coords[nbr[j]])).operator -((coords[nbr[j+1]]).operator *(2))).operator *(0.16));
                        nik=(((coords[i].operator +(coords[nbr[j+1]])).operator -((coords[nbr[j]]).operator *(2))).operator *(0.16)).cross(n);
                        value=(nji.operator *(0.416)).operator +(nik.operator *(0.16));

                        value=value.operator -(normals[i].operator *(value.dot(normals[i])));
                        pesi.push_back(value);
                        area+=m.poly_area(facce[0]);

                        
                    }
                    else if(j==nbr.size()-1)
                    {
                        vec3d n=normals[facce[0]];
                        nij=(((coords[i].operator +(coords[nbr[j]])).operator -((coords[nbr[j-1]]).operator *(2))).operator *(0.16)).cross(n);
                        nhi=n.cross(((coords[i].operator +(coords[nbr[j-1]])).operator -((coords[nbr[j]]).operator *(2))).operator *(0.16));
                        value=(nij.operator *(0.416)).operator +(nhi.operator *(0.16));
                        value=value.operator -(normals[i].operator *(value.dot(normals[i])));
                        pesi.push_back(value);
                        for(int h=0;h<pesi.size();++h)
                        {
                            row=i*3;
                            
                            entries.push_back(Entry(row,nbr[h],pesi.at(h)[0]/area));
                            entries.push_back(Entry(row+1,nbr[h],pesi.at(h)[1]/area));
                            entries.push_back(Entry(row+2,nbr[h],pesi.at(h)[2]/area));
                        }
                        pesi.clear();
                    }else
                    {
                        if(m.poly_contains_vert(facce[1],nbr[j+1]))
                        {
                            Njh=normals[facce[0]];
                            Njk=normals[facce[1]];
                        }else
                        {
                            Njh=normals[facce[1]];
                            Njk=normals[facce[0]];
                        }
                        nij=(((coords[i].operator +(coords[nbr[j]])).operator -((coords[nbr[j-1]]).operator *(2))).operator *(0.16)).cross(Njh);
                        nji=Njk.cross(((coords[i].operator +(coords[nbr[j]])).operator -((coords[nbr[j+1]]).operator *(2))).operator *(0.16));
                        nik=(((coords[i].operator +(coords[nbr[j+1]])).operator -((coords[nbr[j]]).operator *(2))).operator *(1/6)).cross(Njk);
                        nhi=Njh.cross(((coords[i].operator +(coords[nbr[j-1]])).operator -((coords[nbr[j]]).operator *(2))).operator *(0.16));
                        value=(((nij.operator +(nji)).operator *(0.416)).operator +(nik.operator *(0.16))).operator +(nhi.operator *(0.16));
                        value=value.operator -(normals[i].operator *(value.dot(normals[i])));
                        pesi.push_back(value);
                        area+=m.poly_area(facce[1]);
                    }
                }
            */}else{
                nbr=m.vert_ordered_vert_ring(i);
                std::vector<vec3d> pesi(nbr.size());
                vertex_contribute=vec3d(0,0,0);

                area=0.0;
                for (uint j=0;j<nbr.size();++j)
                {
                    int h=(j==0)? nbr.size()-1 : j-1;
                    int k=(j==nbr.size()-1)? 0:j+1;
                    sx={nbr[h],nbr[j],i};
                    dx={nbr[k],nbr[j],i};
                    uint Tjh=m.poly_id(sx);
                    uint Tjk=m.poly_id(dx);
                    area+=m.poly_area(Tjk)/3;

                    Njh=normals[Tjh];
                    Njk=normals[Tjk];


                    first=coords[i].operator *(0.16);
                    second=coords[nbr[j]].operator *(0.16);
                    third=coords[nbr[h]].operator *(0.33);
                    nij=first.operator +(second);
                    nij=nij.operator -(third);
                    nij=nij.cross(Njh);




                    third=coords[nbr[k]].operator *(0.33);
                    nji=first.operator +(second);
                    nji=nji.operator -(third);
                    nji=Njk.cross(nji);


                    second=coords[nbr[k]].operator *(0.16);
                    third=coords[nbr[j]].operator *(0.33);
                    nik=first.operator +(second);
                    nik=nik.operator -(third);
                    nik=nik.cross(Njk);



                    second=coords[nbr[h]].operator *(0.16);
                    nhi=first.operator +(second);
                    nhi=nhi.operator -(third);
                    nhi=Njh.cross(nhi);


                    value=nij.operator +(nji);
                    value=value.operator *(0.416);



                    nik=nik.operator *(0.16);
                    nhi=nhi.operator *(0.16);
                    value=value.operator +(nik);
                    value=value.operator +(nhi);

                    value=vertex_normals[i].cross(value);
                    value=value.cross(vertex_normals[i]);






                    pesi[j]=value;
                    if(j==nbr.size()-1)
                    {
                        for(uint id=0;id<pesi.size();++id)
                        {
                            row=i*3;
                            entries.push_back(Entry(row,nbr[id],pesi.at(id)[0]/area));
                            entries.push_back(Entry(row+1,nbr[id],pesi.at(id)[1]/area));
                            entries.push_back(Entry(row+2,nbr[id],pesi.at(id)[2]/area));
                        }
                    }
                }
            }
        }
        
        G.setFromTriplets(entries.begin(), entries.end());

        return G;
    }
    else if(mode==2)
    {
        Eigen::SparseMatrix<double> G(m.num_verts()*3, m.num_verts());
        std::vector<Entry> entries;
        uint Nf=m.num_polys();
        uint Nv=m.num_verts();
        std::vector<uint> nbr;
        std::vector<uint> sx;
        std::vector<uint> dx;
        vec3d nij;
        vec3d nji;
        vec3d nik;
        vec3d nhi;
        vec3d Njh;
        vec3d Njk;
        vec3d value;
        double lij;
        double lji;
        double lik;
        double lhi;

        vec3d vertex_contribute;
        vec3d first;
        vec3d second;
        vec3d third;
        double area=0.0;
        std::vector<vec3d> normals(Nf);
        std::vector<vec3d> vertex_normals(Nv);
        std::vector<vec3d> coords=m.vector_verts();
        int row=0;

        for(uint i=0; i<Nf;++i)
        {
            normals[i]=m.poly_data(i).normal;
            normals[i].normalize();
        }

        for(uint j=0;j<Nv;++j)
        {
            vertex_normals[j]=m.vert_data(j).normal;
            vertex_normals[j].normalize();
        }

        for (uint i=0;i<Nv;++i)
        {
            /*if(m.vert_is_boundary(i)==1)
            {
                nbr=m.vert_ordered_vert_ring(i);
                area=0.0;
                for (int j=0;j<nbr.size();++j)
                {
                    int eid=m.edge_id(i,nbr[j]);
                    std::vector<uint>facce=m.adj_e2p(eid);
                    if(j==0)
                    {
                        vec3d n=normals[facce[0]];
                        nji=n.cross(((coords[i].operator +(coords[nbr[j]])).operator -((coords[nbr[j+1]]).operator *(2))).operator *(0.16));
                        nik=(((coords[i].operator +(coords[nbr[j+1]])).operator -((coords[nbr[j]]).operator *(2))).operator *(0.16)).cross(n);
                        value=(nji.operator *(0.416)).operator +(nik.operator *(0.16));

                        value=value.operator -(normals[i].operator *(value.dot(normals[i])));
                        pesi.push_back(value);
                        area+=m.poly_area(facce[0]);


                    }
                    else if(j==nbr.size()-1)
                    {
                        vec3d n=normals[facce[0]];
                        nij=(((coords[i].operator +(coords[nbr[j]])).operator -((coords[nbr[j-1]]).operator *(2))).operator *(0.16)).cross(n);
                        nhi=n.cross(((coords[i].operator +(coords[nbr[j-1]])).operator -((coords[nbr[j]]).operator *(2))).operator *(0.16));
                        value=(nij.operator *(0.416)).operator +(nhi.operator *(0.16));
                        value=value.operator -(normals[i].operator *(value.dot(normals[i])));
                        pesi.push_back(value);
                        for(int h=0;h<pesi.size();++h)
                        {
                            row=i*3;

                            entries.push_back(Entry(row,nbr[h],pesi.at(h)[0]/area));
                            entries.push_back(Entry(row+1,nbr[h],pesi.at(h)[1]/area));
                            entries.push_back(Entry(row+2,nbr[h],pesi.at(h)[2]/area));
                        }
                        pesi.clear();
                    }else
                    {
                        if(m.poly_contains_vert(facce[1],nbr[j+1]))
                        {
                            Njh=normals[facce[0]];
                            Njk=normals[facce[1]];
                        }else
                        {
                            Njh=normals[facce[1]];
                            Njk=normals[facce[0]];
                        }
                        nij=(((coords[i].operator +(coords[nbr[j]])).operator -((coords[nbr[j-1]]).operator *(2))).operator *(0.16)).cross(Njh);
                        nji=Njk.cross(((coords[i].operator +(coords[nbr[j]])).operator -((coords[nbr[j+1]]).operator *(2))).operator *(0.16));
                        nik=(((coords[i].operator +(coords[nbr[j+1]])).operator -((coords[nbr[j]]).operator *(2))).operator *(1/6)).cross(Njk);
                        nhi=Njh.cross(((coords[i].operator +(coords[nbr[j-1]])).operator -((coords[nbr[j]]).operator *(2))).operator *(0.16));
                        value=(((nij.operator +(nji)).operator *(0.416)).operator +(nik.operator *(0.16))).operator +(nhi.operator *(0.16));
                        value=value.operator -(normals[i].operator *(value.dot(normals[i])));
                        pesi.push_back(value);
                        area+=m.poly_area(facce[1]);
                    }
                }
            }else{*/

            nbr=m.vert_ordered_vert_ring(i);
            std::vector<vec3d> pesi(nbr.size());
            vertex_contribute=vec3d(0,0,0);

            area=0.0;
            for(uint j=0;j<nbr.size();++j)
            {
                int h=(j==0)? nbr.size()-1 : j-1;
                int k=(j==nbr.size()-1)? 0:j+1;
                sx={nbr[h],nbr[j],i};
                dx={nbr[k],nbr[j],i};
                uint Tjh=m.poly_id(sx);
                uint Tjk=m.poly_id(dx);
                area+=m.poly_area(Tjk)/3;

                Njh=normals[Tjh];
                Njk=normals[Tjk];

                first=coords[i].operator *(0.16);
                second=coords[nbr[j]].operator *(0.16);
                third=coords[nbr[h]].operator *(0.33);
                nij=first.operator +(second);
                nij=nij.operator -(third);
                //third=nij;
                lij=nij.length();
                nij=nij.operator /(2);
                nij=nij.cross(Njh);
                //nij=nij.operator *(lij);
                nij=vertex_normals[i].cross(nij);
                nij=nij.cross(vertex_normals[i]);
                //nij=nij.cross(third);

                third=coords[nbr[k]].operator *(0.33);
                nji=first.operator +(second);
                nji=nji.operator -(third);
                lji=nji.length();
                //third=nji;
                //third=third.operator -();
                nji=nji.operator /(2);
                nji=Njk.cross(nji);
                //nji=nji.operator *(lji);
                nji=vertex_normals[i].cross(nji);
                nji=nji.cross(vertex_normals[i]);
                //nji=nji.cross(third);

                second=coords[nbr[k]].operator *(0.16);
                third=coords[nbr[j]].operator *(0.33);
                nik=first.operator +(second);
                nik=nik.operator -(third);
                lik=nik.length();
                //third=nik;
                nik=nik.operator /(2);
                nik=nik.cross(Njk);
                //nik=nik.operator *(lik);
                nik=vertex_normals[i].cross(nik);
                nik=nik.cross(vertex_normals[i]);
                //nik=nik.cross(third);

                second=coords[nbr[h]].operator *(0.16);
                nhi=first.operator +(second);
                nhi=nhi.operator -(third);
                lhi=nhi.length();
                //third=nhi;
                //third=third.operator -();
                nhi=nhi.operator /(2);
                nhi=Njh.cross(nhi);
                //nhi=nhi.operator *(lhi);
                //nhi=nhi.cross(third);
                nhi=vertex_normals[i].cross(nhi);
                nhi=nhi.cross(vertex_normals[i]);

                value=nij.operator +(nji);
                value=value.operator *(0.416);

                nik=nik.operator *(0.16);
                nhi=nhi.operator *(0.16);
                value=value.operator +(nik);
                value=value.operator +(nhi);

                /*value=vertex_normals[i].cross(value);
                value=value.cross(vertex_normals[i]);*/

                /*provv=value.dot(vertex_normals[i]);
                covalue=vertex_normals[i].operator *(provv);
                value=value-covalue;*/

                pesi[j]=value;
                if(j==nbr.size()-1)
                {
                    for(uint id=0;id<pesi.size();++id)
                    {
                        row=i*3;
                        entries.push_back(Entry(row,nbr[id],pesi.at(id)[0]/area));
                        entries.push_back(Entry(row+1,nbr[id],pesi.at(id)[1]/area));
                        entries.push_back(Entry(row+2,nbr[id],pesi.at(id)[2]/area));

                    }
                /*vertex_contribute=vertex_normals[i].cross(vertex_contribute);
                vertex_contribute=vertex_contribute.cross(vertex_normals[i]);*/


                /*entries.push_back(Entry(row,i,vertex_contribute[0]/area));
                entries.push_back(Entry(row+1,i,vertex_contribute[1]/area));
                entries.push_back(Entry(row+2,i,vertex_contribute[2]/area));*/

                }
            }
        }

        G.setFromTriplets(entries.begin(), entries.end());

        return G;
    }
}



