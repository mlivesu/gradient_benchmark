#include "computations.h"
#include "picking.h"
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <cinolib/symbols.h>
#include <cinolib/color.h>
#include <cinolib/geometry/triangle.h>
#include <iostream>
#include <cinolib/gradient.h>
#include <cinolib/point_inside_mesh.h>

extern COMPUTATIONS computations;

using namespace cinolib;
using namespace Eigen;
typedef Eigen::Triplet<double> Entry;
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
bool poly_has_vert_on_srf (const DrawableTrimesh<> &m, uint pid)
{
    for(uint vid : m.adj_p2v(pid))
    {
        if(m.vert_is_boundary(vid)) return true;
    }
    return false;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void num_verts_on_boundary(const DrawableTrimesh<> & m, std::map<uint,int> &boundaries)
{
    int count=0;
    for(uint i=0;i<m.num_verts();++i)
    {
        if(m.vert_is_boundary(i))
        {

            boundaries[i]=count;
            ++count;

        }
    }



}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

double relative_error(const vec3d v1, const vec3d v2,const int mode)
{
    double norm1;
    double norm2;
    double angle;
    double diff_norm;
    double err;
    vec3d diff;
    double treshold=1e-8;
    switch(mode)
    {
    case 0:

    {
        diff=v1.operator -(v2);
        norm1=diff.length();
        norm2=(v1.length()+v2.length())/2;
        err=(norm2>treshold)? norm1/norm2 : treshold;

        break;
    }
    case 1:

    {
        norm1=v1.length();
        norm2=v2.length();
        diff_norm=abs(norm1-norm2);
        double tmp=(norm1 +norm2)/2;
        err=(tmp>treshold)? diff_norm/tmp  : treshold;


        break;
    }
    case 2:

    {

        if(v1.length()<=1e-16 || v2.length()<=1e-16 )
        {
            err=1e-16;
        }else
        {
            err=pow(v1.angle_rad(v2),2);
        }


        break;
    }


    }

    return err;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double absolute_error(const vec3d v1, const vec3d v2, const int mode)
{
    double norm1;
    double norm2;
    double err;
    vec3d diff;
    switch(mode)
    {
    case 0:

        diff=v1.operator -(v2);
        err=diff.length_squared();
        break;
    case 1:

        norm1=v1.length();
        norm2=v2.length();
        err=pow(norm1-norm2,2);


        break;
    case 2:

        if(v1.length()<=1e-16 || v2.length()<=1e-16 )
        {
            err=1e-16;
        }else
        {
            err=pow(v1.angle_rad(v2),2);
        }



        break;
    }

    return err;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double scale_function(double x, double md, double Md)
{


    double amplitude=Md-md;
//    double delta=amplitude/25;
//    double Mdtmp=md+scale_factor*delta;

    x=(x-md)/(Md-md);

//    if (x>=1){x=0.999;}
//    if (x<=0){x=1e-10;}
    return x;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double find_max_norm(const DrawableVectorField &V)
{
    int N=V.rows()/3;
    double max=0;
    for(int i=0;i<N;++i)
    {
        max=std::max(max,V.vec_at(i).length());
    }
    return max;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void find_max_min_values (const ScalarField f, double &max, double &min)
{
    int N=f.rows();
    max=0;
    min=0;

    for (int i=0;i<N;++i)
    {
        max=std::max(max,f[i]);
        min=std::min(min,f[i]);

    }

}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void pid_contribute(const std::vector<imaginary_vertex> &nbr, std::map<uint, int> boundaries, std::vector<std::pair<uint, vec3d> > &pesi, int &offset, vec3d &n)
{
    int entry;
    for(uint off=0; off<3; ++off)
    {

        uint  curr = nbr[(off+1)%3].real;


        vec3d u    = nbr[(off+2)%3].p - nbr[(off+1)%3].p;
        vec3d v    = nbr[(off+1)%3].p - nbr[off].p;
        vec3d u_90 = u.cross(n); u_90.normalize();
        vec3d v_90 = v.cross(n); v_90.normalize();

        if(curr>offset)
        {
            curr-=offset;

            entry=boundaries[curr];
            pesi.push_back(std::make_pair(offset+entry,u_90 + v_90));
        }else
        {
            pesi.push_back(std::make_pair(curr,u_90 + v_90));
        }


    }

}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
vec3d generic_contribute(const vec3d vi,const vec3d vj,const vec3d vh, const vec3d vk, const vec3d Njh, const vec3d Njk,int mode)

{
    vec3d vij=vi.operator -(vj);
    vec3d vji=vij.operator -();

    vec3d vjh=vj.operator -(vh);
    vec3d vhj=vjh.operator -();

    vec3d vkj=vk.operator -(vj);
    vec3d vjk=vkj.operator -();

    vec3d contribute;

    vec3d vij_90=vij.cross(Njh);
    vec3d vjh_90=vjh.cross(Njh);
    vec3d vji_90=vji.cross(Njk);
    vec3d vkj_90=vkj.cross(Njk);

    switch(mode)
    {
    case 0:
    {

        contribute = vkj_90 + vji_90 + vij_90 + vjh_90;

    }
        break;
    case 1:
    {

        contribute = vkj_90 + vji_90;


    }
        break;
    case 2 :
    {

        contribute = vij_90 + vjh_90;

    }
        break;
    }

    return contribute;

}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
vec3d vi_contribute(const vec3d vi,const vec3d vj,const vec3d vh, const vec3d vk, const vec3d Njh, const vec3d Njk,int mode)
{
    vec3d vij=vi.operator -(vj);
    vec3d vih=vi.operator -(vh);
    vec3d vhi=vh.operator -(vi);
    vec3d vji=vj.operator -(vi);
    vec3d vik=vi.operator -(vk);
    vec3d vki=vik.operator -();


    vec3d vij_90=vij.cross(Njh);
    vec3d vhi_90=vhi.cross(Njh);
    vec3d vji_90=vji.cross(Njk);
    vec3d vik_90=vik.cross(Njk);
    vec3d vki_90=vki.cross(Njk);
    vec3d contribute;

    switch(mode)
    {

    case 0:
    {

    }
        break;
    case 1:
    { contribute= vji_90 + vik_90;
    }
        break;
    case 2:
    {

        contribute= vhi_90 + vij_90;


    }
        break;

    }

    return contribute;



}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
std::vector<imaginary_vertex> nbr_for_boundaries(const DrawableTrimesh<> &m, uint vid)
{
    std::vector<imaginary_vertex> nbr;
    int N=m.num_verts();
    std::vector<uint> tmp=m.vert_ordered_vert_ring(vid);
    double l=m.edge_avg_length();
    vec3d prev=m.vert(tmp[0]);
    vec3d next=m.vert(tmp[tmp.size()-1]);
    vec3d dir_prev;
    vec3d dir_curr;
    vec3d dir_next;
    if(prev[0]==next[0]||prev[1]==next[1])
    {
        dir_prev=-(m.vert(vid)-prev).cross(m.vert_data(tmp[0]).normal);
        dir_curr=(prev-m.vert(vid)).cross(m.vert_data(vid).normal);
        dir_next=(m.vert(vid)-next).cross(m.vert_data(tmp[tmp.size()-1]).normal);
        dir_prev.normalize();
        dir_prev*=l;
        dir_prev+=prev;
        dir_curr.normalize();
        dir_curr*=l;
        dir_prev+=m.vert(vid);
        dir_next.normalize();
        dir_next*=l;
        dir_next+=next;
    }else
    {
        dir_prev=(m.vert(vid)-prev).cross(-m.vert_data(tmp[0]).normal);
        dir_next=(m.vert(vid)-next).cross(m.vert_data(tmp[tmp.size()-1]).normal);
        dir_curr=dir_prev+dir_next;
        dir_curr/=2;
        dir_prev.normalize();
        dir_prev*=l;
        dir_prev+=prev;
        dir_curr.normalize();
        dir_curr*=l;
        dir_prev+=m.vert(vid);
        dir_next.normalize();
        dir_next*=l;
        dir_next+=next;
    }

    imaginary_vertex v0;
    imaginary_vertex v1;
    imaginary_vertex v2;

    v0.p=dir_prev;
    v0.real=N+tmp[0];
    nbr.push_back(v0);
    v1.p=dir_curr;
    v1.real=N+vid;
    nbr.push_back(v1);
    v2.p=dir_next;
    v2.real=N+tmp[tmp.size()-1];
    nbr.push_back(v2);


    return nbr;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

Eigen::SparseMatrix<double> GG_gradient_matrix(const DrawableTrimesh<> & m,int mode)
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
    else if(mode==1)
    {
        std::map<uint,int> boundaries;
        num_verts_on_boundary(m,boundaries);
        int B=boundaries.size();
        int N=m.num_verts();
        int size=N+B;

        Eigen::SparseMatrix<double> G(m.num_verts()*3, size);
        std::vector<Entry> entries;
        std::vector<std::pair<uint,vec3d>> pesi;

        for (uint i=0;i<N;++i)
        {
            double area=0.0;
            pesi.clear();

            for(uint pid : m.adj_v2p(i))
            {
                area+=2*m.poly_area(pid);
                vec3d n     = m.poly_data(pid).normal;

                for(uint off=0; off<m.verts_per_poly(pid); ++off)
                {
                    uint  prev = m.poly_vert_id(pid,off);
                    uint  curr = m.poly_vert_id(pid,(off+1)%m.verts_per_poly(pid));
                    uint  next = m.poly_vert_id(pid,(off+2)%m.verts_per_poly(pid));
                    vec3d u    = m.vert(next) - m.vert(curr);
                    vec3d v    = m.vert(curr) - m.vert(prev);
                    vec3d u_90 = u.cross(n);
                    vec3d v_90 = v.cross(n);

                    pesi.push_back(std::make_pair(curr,u_90 + v_90));

                }
            }

            if(m.vert_is_boundary(i))
            {
                std::vector<imaginary_vertex>nbr_aux=nbr_for_boundaries(m,i);
                std::vector<uint> nbr=m.vert_ordered_vert_ring(i);
                imaginary_vertex v_0;
                imaginary_vertex v_n;
                imaginary_vertex v_i;
                v_i.p=m.vert(i);
                v_i.real=i;
                v_0.p=m.vert(nbr[0]);
                v_0.real=nbr[0];
                v_n.p=m.vert(nbr[nbr.size()-1]);
                v_n.real=nbr[nbr.size()-1];

                for(uint j=0;j<4;++j)
                {
                    imaginary_vertex v0=(j==0)?v_0:nbr_aux[j-1];
                    imaginary_vertex v2=(j==3)?v_n:nbr_aux[j+1];

                    vec3d tmp0 =v0.p - v_i.p;
                    vec3d tmp1=v2.p-v_i.p;
                    vec3d n=tmp0.cross(tmp1);
                    area+=n.length();
                    std::vector<imaginary_vertex>triangle;
                    triangle.push_back(v0);
                    triangle.push_back(v2);
                    triangle.push_back(v_i);
                    n.normalize();

                    pid_contribute(triangle,boundaries,pesi,N,n);

                }
            }

            int  row=i*3;

            for(int id=0;id<pesi.size();++id)
            {

                entries.push_back(Entry(row,pesi[id].first,pesi[id].second[0]/area));
                entries.push_back(Entry(row+1,pesi[id].first,pesi[id].second[1]/area));
                entries.push_back(Entry(row+2,pesi[id].first,pesi[id].second[2]/area));
            }
        }


        G.setFromTriplets(entries.begin(), entries.end());

        return G;
    }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
DrawableVectorField from_f2v(DrawableVectorField & W, const DrawableTrimesh<> & m)

{

    int Nv=m.num_verts();
    std::vector<vec3d> grad_v(Nv);
    vec3d provv=vec3d(0,0,0);
    vec3d avg=vec3d(0,0,0);

    std::vector<uint> nbr;
    DrawableVectorField Wv;
    Wv=DrawableVectorField(m,false);


    for (int i=0;i<Nv;++i)
    {
        nbr=m.adj_v2p(i);
        provv=vec3d(0,0,0);
        avg=vec3d(0,0,0);
        vec3d pos=m.vert(i);
        double wgt=0;
        for(int j=0;j<nbr.size();++j)
        {
            uint pid=nbr[j];
            double weight=1/(m.poly_centroid(pid)-pos).length();
            provv=W.vec_at(pid)*weight;
            wgt+=weight;
            avg+=provv;

        }
        avg/=wgt;
        Wv.set(i,avg);
    }
    return Wv;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
DrawableVectorField  compute_ground_truth (const DrawableTrimesh<> &m, const double a, const double b, const double c, const int mode, const int method)
{
    DrawableVectorField V;
    int N=0;
    std::vector<vec3d> coords;
    double A=a*10;
    double B=b/100;
    if(method==0)
    {
        V=DrawableVectorField(m,true);
        N=m.num_polys();
        for (int j=0;j<N;++j)
        {
            coords.push_back(m.poly_centroid(j));
        }
    }else
    {

        V=DrawableVectorField(m,false);
        N=m.num_verts();
        coords=m.vector_verts();
    }


    cinolib::Color col;
    col=col.BLACK();
    switch(mode)
    {

    case 0:
        for (int i=0;i<N;++i)
        {
            vec3d pos=coords[i];
            vec3d grad=vec3d((A*sin(3*b*pos[0])+A*3*b*pos[0]*cos(3*b*pos[0]))*sin(3*b*pos[1]*pos[1]),
                    6*b*A*pos[0]*pos[1]*sin(3*b*pos[0])*cos(3*b*pos[1]*pos[1]),0);
            V.set(i,grad);
        }
        break;
    case 1:
        for (int i=0;i<N;++i)
        {
            vec3d pos=coords[i];
            vec3d grad=vec3d(2*b*(pos[0]-0.5),2*c*(pos[1]-0.5),0);
            V.set(i,grad);
        }
        break;
    case 2:
        for(int i=0;i<N;++i)
        {
            vec3d pos=coords[i];
            vec3d grad=vec3d(A*B*cos(B*pos[0])*cos(B*pos[1]),-A*B*sin(B*pos[0])*sin(B*pos[1]),0);
            V.set(i,grad);
        }
        break;
    case 3:
        for (int i=0;i<N;++i)
        {

            vec3d grad=vec3d(1,1,0);
            V.set(i,grad);
        }
        break;
    }


    V.set_arrow_color(col);
    return V;

}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
ScalarField get_scalar_field(const DrawableTrimesh<> & m, const double a, const double b, const double c, const int mode)
{
    std::vector<vec3d> v=m.vector_verts();
    int Nv= m.num_verts();
    double A=a*10;
    double B=b/100;
    std::vector<double> data(Nv);
    switch (mode)
    {
    case 0:

        for(int i=0;i<Nv;++i)
        {
            vec3d w=v[i];
            double val=A*w[0]*sin(3*b*w[0])*sin(3*b*w[1]*w[1]);
            data[i]=val;
        }
        break;
    case 1:
        for(int i=0;i<Nv;++i)
        {
            vec3d w=v[i];
            double val=b*pow(w[0]-0.5,2)+c*pow(w[1]-0.5,2);
            data[i]=val;
        }
        break;
    case 2:
        for(int i=0;i<Nv;++i)
        {
            vec3d w=v[i];
            double val=A*sin(B*w[0])*cos(B*w[1]);
            data[i]=val;
        }
        break;
    case 3:

        for(int i=0;i<Nv;++i)
        {
            vec3d w=v[i];
            double val=w[0]+w[1];
            data[i]=val;
        }

        break;
    }


    return ScalarField(data);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
ScalarField scalar_field_with_boundaries(const DrawableTrimesh<> &m, const double a, const double b, const double c, const int mode, std::map<uint, int> &boundaries)
{
    std::vector<vec3d> v=m.vector_verts();
    int Nv= m.num_verts();
    num_verts_on_boundary(m,boundaries);
    int B=boundaries.size();
    int s=Nv+B;
    double A=a/10;
    std::vector<double> data(s);
    switch (mode)
    {
    case 0:

        for(int i=0;i<Nv;++i)
        {
            vec3d w=v[i];
            double val=A*w[0]*sin(3*b*w[0])*sin(3*b*w[1]*w[1]);
            data[i]=val;
        }
        break;
    case 1:
        for(int i=0;i<Nv;++i)
        {
            vec3d w=v[i];
            double val=b*pow(w[0]-0.5,2)+c*pow(w[1]-0.5,2);
            data[i]=val;
        }
        break;
    case 2:
        for(int i=0;i<Nv;++i)
        {
            vec3d w=v[i];
            double val=A*sin(b*w[0])*cos(b*w[1]);
            data[i]=val;
            if(m.vert_is_boundary(i))
            {
                int entry=boundaries[i];
                data[Nv+entry]=val;
            }

        }
        break;
    case 3:

        for(int i=0;i<Nv;++i)
        {
            vec3d w=v[i];
            double val=w[0]+w[1];
            data[i]=val;
        }

        break;
    }


    return ScalarField(data);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

DrawableVectorField compute_field(DrawableTrimesh<> & m, ScalarField & f, const int mode, bool handle_boundary)
{

    DrawableVectorField V;
    DrawableVectorField W;
    Eigen::SparseMatrix<double> G;
    cinolib::Color c;
    c=c.RED();
    switch(mode)
    {
    case 0:
    {
        G=GG_gradient_matrix(m,0);
        V=DrawableVectorField(m,true);
        V=G*f;
        V.set_arrow_color(c);
    }

        break;
    case 1:
    {
        if(handle_boundary)
            G=GG_gradient_matrix(m,1);
        else
            G=gradient_matrix(m,0);

        V=DrawableVectorField(m,false);
        V=G*f;
        V.set_arrow_color(c);
    }
        break;
    case 2:
    {

        V=compute_FEM(m,f,handle_boundary);
        V.set_arrow_color(c);
    }
        break;
    case 3:
    {
        V=compute_quadratic_regression(m,f,handle_boundary);
        V.set_arrow_color(c);

    }
        break;

    case 4:
    {
        V=compute_quadratic_regression_centroids(m,f,handle_boundary);
        V.set_arrow_color(c);
    }
        break;

    }
    return V;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
ScalarField estimate_error(const DrawableVectorField & GT, const DrawableVectorField & V, const DrawableTrimesh<> &m, int mode, int method, int relative, const int type_of_vertices)
{
    std::vector<double> data;
    double err=0;
    double max=0;
    double min=1e+10;
    double avg=0;
    int count=0;
    if(method!=0)
    {
        int N=m.num_verts();
        switch(type_of_vertices)
        {
        case 0:
        {
            for(int i=0;i<N;++i)
            {
                err=(relative==0)?absolute_error(GT.vec_at(i),V.vec_at(i),mode):relative_error(GT.vec_at(i),V.vec_at(i),mode);
                data.push_back(err);
                avg+=err;
                max=std::max(err,max);
                min=std::min(err,min);
            }
            avg=avg/N;
            if(mode==2)
            {
                avg=sqrt(avg);
            }else
            {
                avg/=2;
            }
            std::cout<<"Err MAX="<<max<<std::endl;
            std::cout<<"Err MIN="<<min<<std::endl;
            std::cout<<"MEAN Err="<<avg<<std::endl;

        } break;
        case 1:
        {
            for(int i=0;i<N;++i)
            {
                if(m.vert_is_boundary(i)||m.vert_data(i).marked)
                    data.push_back(1e-10);
                else
                {
                    err=(relative==0)?absolute_error(GT.vec_at(i),V.vec_at(i),mode):relative_error(GT.vec_at(i),V.vec_at(i),mode);

                        avg+=err;
                        max=std::max(err,max);
                        min=std::min(err,min);
                        ++count;
                        data.push_back(err);

                }

            }
            avg=avg/count;
            if(mode==2)
            {
                avg=sqrt(avg);
            }else
            {
                avg/=2;
            }
            std::cout<<"Err MAX="<<max<<std::endl;
            std::cout<<"Err MIN="<<min<<std::endl;
            std::cout<<"MEAN Err="<<avg<<std::endl;

        } break;
        case 2:
        {
            for(int i=0;i<N;++i)
            {
                if(m.vert_is_boundary(i))
                {
                     ++count;

                    err=(relative==0)?absolute_error(GT.vec_at(i),V.vec_at(i),mode):relative_error(GT.vec_at(i),V.vec_at(i),mode);

                    data.push_back(err);
                    avg+=err;
                    max=std::max(err,max);
                    min=std::min(err,min);
                }
                else
                {

                    data.push_back(1e-10);
                }

            }
            avg=avg/count;
            if(mode==2)
            {
                avg=sqrt(avg);
            }else
            {
                avg/=2;
            }
            std::cout<<"Err MAX="<<max<<std::endl;
            std::cout<<"Err MIN="<<min<<std::endl;
            std::cout<<"MEAN Err="<<avg<<std::endl;

        }break;

        }


    }
    else
    {
        data=dual_error(m,GT,V,mode,relative);
    }


    return ScalarField(data);

}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
std::vector<double> dual_error(const DrawableTrimesh<> &m, const DrawableVectorField &GT, const DrawableVectorField &V, const int mode,const int relative)
{



    std::vector<vec3d>             dual_verts;
    std::vector<std::vector<uint>> dual_faces;
    dual_mesh(m, dual_verts,dual_faces,true);
    computations.dual_m=DrawablePolygonmesh<>(dual_verts,dual_faces);

    std::vector<double> data;
    double err=0;
    double max=0;
    double min=1e+5;
    double avg=0;
    double count=0;
    int M=m.num_polys();

    switch(relative)
    {
    case 0:
        for(int i=0;i<computations.dual_m.num_verts();++i)
        {
            if(i<M)
            {

                err=absolute_error(GT.vec_at(i),V.vec_at(i),mode);
                data.push_back(err);
                avg+=err;
                ++count;
                max=std::max(err,max);
                min=std::min(err,min);

            }else
            {
                data.push_back(0);
            }



        }
        avg=avg/count;
        std::cout<<"Err MAX="<<max<<std::endl;
        std::cout<<"Err MIN="<<min<<std::endl;
        std::cout<<"MEAN Err="<<avg<<std::endl;
        break;
    case 1:
        for(int i=0;i<computations.dual_m.num_verts();++i)
        {
            if(i<M)
            {

                err=relative_error(GT.vec_at(i),V.vec_at(i),mode);
                data.push_back(err);
                avg+=err;
                ++count;
                max=std::max(err,max);
                min=std::min(err,min);
            }else
            {
                data.push_back(0);
            }


        }
        avg=avg/count;
        std::cout<<"Err MAX="<<max<<std::endl;
        std::cout<<"Err MIN="<<min<<std::endl;
        std::cout<<"MEAN Err="<<avg<<std::endl;

        break;
    }



    return data;




}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
DrawableVectorField compute_quadratic_regression(DrawableTrimesh<> &m, const ScalarField & f, bool handle_boundaries)
{
    int nv=m.num_verts();
    std::vector<int> rank (nv);
    std::vector<uint> nbr;
    vec3d vert;
    Eigen::VectorXd X(6);
    DrawableVectorField V=DrawableVectorField(m,false);
    std::vector<imaginary_vertex> nbr_aux;

    for (int i=0;i<nv;++i)
    {
        nbr=m.vert_ordered_vert_ring(i);
        nbr.push_back(i);
        vert=m.vert(i);
        if(nbr.size()<6)
        {
            std::set<uint> two_ring=m.vert_n_ring(i,2);
            nbr.clear();
            for(uint k : two_ring)
            {
                nbr.push_back(k);

            }
            nbr.push_back(i);
        }

        if(m.vert_is_boundary(i) && handle_boundaries)
        {
            nbr_aux=nbr_for_boundaries(m,i);
        }

        Eigen::MatrixXd coeff(nbr.size()+nbr_aux.size(),6);
        Eigen::VectorXd b(nbr.size()+nbr_aux.size());
        int s=nbr.size();
        for (int j=0;j<nbr.size()+nbr_aux.size();++j)
        {
            vec3d pos;
            double wgt;
            if(j<nbr.size())
            {
                pos=m.vert(nbr[j]);
                wgt=(nbr[j]==i)? 1:1/(pos-vert).length_squared();
                b(j)=f[nbr[j]]*wgt;
            }else
            {
                pos=nbr_aux[j-s].p;
                wgt=1/(pos-vert).length_squared();
                uint curr=nbr_aux[j-s].real-nv;
                b(j)=f[curr]*wgt;
            }

            coeff(j,0)=pow(pos[0],2)*wgt;
            coeff(j,1)=pos[0]*pos[1]*wgt;
            coeff(j,2)=pow(pos[1],2)*wgt;
            coeff(j,3)=pos[0]*wgt;
            coeff(j,4)=pos[1]*wgt;
            coeff(j,5)=1*wgt;
        }
        Eigen::MatrixXd coeffT=Transpose<Eigen::MatrixXd>(coeff);
        Eigen::MatrixXd A=coeffT*coeff;
        Eigen::VectorXd B=coeffT*b;
        Eigen::ColPivHouseholderQR<MatrixXd> dec(A);
        X=dec.solve(B);
        vec3d sol= vec3d (2*X(0)*vert.x()+X(1)*vert.y()+X(3),2*X(2)*vert.y()+X(1)*vert.x()+X(4),0);
        V.set(i,sol);

    }
    computations.rank=rank;
    return V;

}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
DrawableVectorField compute_quadratic_regression_centroids(DrawableTrimesh<> &m, const ScalarField & f, bool handle_boundaries)
{
    int nv=m.num_verts();
    std::vector<int> rank (nv);
    std::vector<uint> nbr;
    std::vector<double> values;
    vec3d vert;
    Eigen::VectorXd X(6);
    DrawableVectorField W=DrawableVectorField(m,false);
    std::vector<imaginary_vertex> nbr_aux;

    for (int i=0;i<nv;++i)
    {

        nbr=m.vert_ordered_vert_ring(i);
        nbr.push_back(i);
        values.clear();
        vert=m.vert(i);
        //        if(m.vert_is_boundary(i) && handle_boundaries)
        //        {
        //            nbr_aux=nbr_for_boundaries(m,i);
        //        }
        int V=nbr.size();
        //        int I=nbr_aux.size();
        //        int size=V+I;


        if(V<6)
        {
            for(uint j=0;j<V;++j)
            {
                for(uint pid : m.adj_v2p(nbr[j]))
                {
                    if(!m.poly_data(pid).marked)
                    {
                        nbr.push_back(pid);
                        std::vector<uint> pid_vertices=m.poly_verts_id(pid);
                        double value=0;
                        for(uint j=0;j<3;++j)
                        {
                            value+=f[pid_vertices[j]];
                        }
                        values.push_back(value/3);
                        m.poly_data(pid).marked=true;

                    }
                }
            }
        }
        int P=values.size();
        int size=P+V;

        Eigen::MatrixXd coeff(size,6);
        Eigen::VectorXd b(size);

        for (int j=0;j<size;++j)
        {
            vec3d pos;
            double wgt;
            if(j<V)
            {
                uint curr=nbr[j];
                pos=m.vert(curr);
                wgt=(curr==i)?:1/(pos-vert).length_squared();
                b(j)=f[curr]*wgt;

            }else
            {
                uint curr=nbr[j];
                pos=m.poly_centroid(curr);
                wgt=1/(pos-vert).length_squared();
                b(j)=values[j-V]*wgt;
            }

            coeff(j,0)=pow(pos[0],2)*wgt;
            coeff(j,1)=pos[0]*pos[1]*wgt;
            coeff(j,2)=pow(pos[1],2)*wgt;
            coeff(j,3)=pos[0]*wgt;
            coeff(j,4)=pos[1]*wgt;
            coeff(j,5)=1*wgt;
        }

        Eigen::MatrixXd coeffT=Transpose<Eigen::MatrixXd>(coeff);
        Eigen::MatrixXd A=coeffT*coeff;
        Eigen::VectorXd B=coeffT*b;
        Eigen::ColPivHouseholderQR<MatrixXd> dec(A);

        X=dec.solve(B);

        vec3d sol= vec3d (2*X(0)*vert.x()+X(1)*vert.y()+X(3),2*X(2)*vert.y()+X(1)*vert.x()+X(4),0);
        W.set(i,sol);


    }
    computations.rank=rank;
    return W;

}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
DrawableVectorField compute_FEM(DrawableTrimesh<> &m, const ScalarField & f, bool handle_bondary)
{


    int Nv=m.num_verts();
    int M=0;
    double delta=0;
    std::vector<uint> nbr;
    vec3d grad;
    std::vector<vec3d> coords=m.vector_verts();
    DrawableVectorField V=DrawableVectorField(m,false);
    Eigen::VectorXd X(2);
    Eigen::VectorXd b(2);
    std::vector<double> values;
    for (int i=0;i<Nv;++i)
    {
        values.clear();
        nbr=m.vert_ordered_vert_ring(i);
        int s=nbr.size();
        if(s<3)
        {
            for(uint j=0;j<s;++j)
            {
                for(uint pid : m.adj_v2p(nbr[j]))
                {
                    if(!m.poly_data(pid).marked)
                    {
                        nbr.push_back(pid);
                        std::vector<uint> pid_vertices=m.poly_verts_id(pid);
                        double value=0;
                        for(uint j=0;j<3;++j)
                        {
                            value+=f[pid_vertices[j]];
                        }
                        values.push_back(value/3);
                        m.poly_data(pid).marked=true;

                    }
                }
            }
        }
        M=s+values.size();
        Eigen::MatrixXd B(M,2);
        Eigen::VectorXd sigma(M);
        for(int j=0;j<M;++j)
        {
            if(j<s)
            {
                vec3d vij=coords[nbr[j]].operator -(coords[i]);
                double wgt=1/vij.length_squared();
                B(j,0)=vij[0]*wgt;
                B(j,1)=vij[1]*wgt;
                sigma(j)=(f[nbr[j]]-f[i])*wgt;
            }
            else
            {
                vec3d vij=m.poly_centroid(nbr[j]).operator -(coords[i]);
                double wgt=1/vij.length_squared();
                B(j,0)=vij[0]*wgt;
                B(j,1)=vij[1]*wgt;
                double value=values[j-s];
                sigma(j)=(value-f[i])*wgt;
            }


            if(j==M-1)
            {
                Eigen::Matrix2d A(2,2);
                Eigen::MatrixXd Bt=Transpose<Eigen::MatrixXd>(B);
                A=Bt*B;
                Eigen::ColPivHouseholderQR<Matrix2d> dec(A);
                b=Bt*sigma;

                X=dec.solve(b);


                grad=vec3d(X(0),X(1),0);
                V.set(i,grad);
            }

        }
    }
    return V;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
DrawableVectorField arrows_normalization(const DrawableTrimesh<> &m, const DrawableVectorField &V, const int mode,const int scale_factor)
{
    int N=0;
    vec3d max=vec3d(0,0,0);
    vec3d min=vec3d(0,0,0);
    DrawableVectorField W;

    double treshold=1e-8;
    if(mode!=0)
    {
        N=m.num_verts();
        W=DrawableVectorField(m,false);


    }else
    {
        N=m.num_polys();
        W=DrawableVectorField(m,true);

    }


    for(int i=0;i<N;++i)
    {
        if(V.vec_at(i).length()>50)
        {
            W.vec_at(i)=V.vec_at(i);
            W.vec_at(i).normalize();
        }else
        {
            W.vec_at(i)=V.vec_at(i);
        }
    }
    return W;

}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
ScalarField heat_map_normalization(const ScalarField &f,double min,double max,double sat_neg,double sat_pos)
{
    ScalarField F=f;
    double upper_bound=sat_pos/1000;
    double lower_bound=sat_neg/1000;
    for(int i=0; i<f.rows(); ++i)
    {
        F[i] =scale_function(F[i],min,max);
        if(F[i]>upper_bound)
            F[i]=0.9999;
        if(F[i]<lower_bound)
            F[i]=1e-16;

    }


    return F;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
vec3d barycentric_coordinates(const vec3d &A,const vec3d &B, const vec3d &C, const vec3d &P)
{
    double area=triangle_area(A,B,C);
    double area_ACP=triangle_area(A,C,P);
    double area_BCP=triangle_area(B,C,P);
    return vec3d(area_BCP/area,area_ACP/area,1-area_BCP/area-area_ACP/area);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void bring_the_field_inside(const DrawableTrimesh<> &m, DrawableTrimesh<> &m_grid,DrawableVectorField &V,DrawableVectorField &W,const int method)
{


    int Nv=m_grid.num_verts();
    W=DrawableVectorField(m_grid,false);
    PointInsideMeshCache<Trimesh<>> cache(m);

    for(uint vid=0;vid<Nv;++vid)
    {

        uint pid;
        vec3d pos=m_grid.vert(vid);
        std::vector<double> wgts;
        cache.locate(m_grid.vert(vid), pid, wgts);

        if(poly_has_vert_on_srf(m,pid))
        {

            m_grid.vert_data(vid).marked=true;
            //continue;
        }//else{
        if(method==0)
        {
            W.set(vid,V.vec_at(pid));
        }else
        {


            std::vector<uint>poly_vertices_id=m.poly_verts_id(pid,true);
            std::vector<vec3d> vertices_coords(poly_vertices_id.size());
            vec3d bary_coords;
            for (uint i=0;i<poly_vertices_id.size();++i)
            {
                vertices_coords[i]= m.vert(poly_vertices_id[i]);
            }

            bary_coords=barycentric_coordinates(vertices_coords[0],vertices_coords[1],vertices_coords[2],pos);

            vec3d interpolated_value=bary_coords[0]*V.vec_at(poly_vertices_id[0])+bary_coords[1]*V.vec_at(poly_vertices_id[1])+bary_coords[2]*V.vec_at(poly_vertices_id[2]);
            W.set(vid,interpolated_value);
        }

        //}


    }

}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
DrawableVectorField compute_PCE(const DrawableTrimesh<> & m, ScalarField &f)
{
    DrawableVectorField W=DrawableVectorField(m,true);

    for(int pid=0; pid<m.num_polys(); ++pid)
    {
        double area = std::max(m.poly_area(pid), 1e-5) * 2.0; // (2 is the average term : two verts for each edge)
        vec3d n     = m.poly_data(pid).normal;
        vec3d contribute(0,0,0);
        for(uint off=0; off<m.verts_per_poly(pid); ++off)
        {
            uint  prev = m.poly_vert_id(pid,off);
            uint  curr = m.poly_vert_id(pid,(off+1)%m.verts_per_poly(pid));
            uint  next = m.poly_vert_id(pid,(off+2)%m.verts_per_poly(pid));
            vec3d u    = m.vert(next) - m.vert(curr);
            vec3d v    = m.vert(curr) - m.vert(prev);
            vec3d u_90 = u.cross(n); u_90.normalize();
            vec3d v_90 = v.cross(n); v_90.normalize();

            vec3d per_vert_sum_over_edge_normals =u_90 * u.length() + v_90 * v.length();

            per_vert_sum_over_edge_normals*=f[curr];

            per_vert_sum_over_edge_normals /= area;
            contribute+=per_vert_sum_over_edge_normals;
        }
        W.set(pid,contribute);

    }
    return W;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
DrawableVectorField compute_AGS(const DrawableTrimesh<> & m, ScalarField &f)
{
    DrawableVectorField W=DrawableVectorField(m,true);
    DrawableVectorField V=DrawableVectorField(m,false);
    W=compute_PCE(m,f);
    V=from_f2v(W,m);

    return V;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Eigen::VectorXd fit_with_quadrics(const DrawableTrimesh<> &m, uint i, ScalarField &f)
{
    std::vector<uint> nbr=m.vert_ordered_vert_ring(i);
    nbr.push_back(i);
    vec3d vert=m.vert(i);
    if(nbr.size()<6)
    {
        std::set<uint> two_ring=m.vert_n_ring(i,2);
        nbr.clear();
        for(uint k : two_ring)
        {
            nbr.push_back(k);

        }
        nbr.push_back(i);
    }
    Eigen::MatrixXd coeff(nbr.size(),6);
    Eigen::VectorXd b(nbr.size());
    int s=nbr.size();
    for (int j=0;j<nbr.size();++j)
    {
        vec3d pos;
        double wgt;
        pos=m.vert(nbr[j]);
        wgt=(nbr[j]==i)? 1:1/(pos-vert).length_squared();
        b(j)=f[nbr[j]]*wgt;

        coeff(j,0)=pow(pos[0],2)*wgt;
        coeff(j,1)=pos[0]*pos[1]*wgt;
        coeff(j,2)=pow(pos[1],2)*wgt;
        coeff(j,3)=pos[0]*wgt;
        coeff(j,4)=pos[1]*wgt;
        coeff(j,5)=1*wgt;
    }
    Eigen::MatrixXd coeffT=Transpose<Eigen::MatrixXd>(coeff);
    Eigen::MatrixXd A=coeffT*coeff;
    Eigen::VectorXd B=coeffT*b;
    Eigen::ColPivHouseholderQR<MatrixXd> dec(A);
    Eigen::VectorXd X=dec.solve(B);

    return X;
}
