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
bool poly_has_vert_on_boundary (const DrawableTrimesh<> &m, uint pid)
{
    for(uint vid : m.adj_p2v(pid))
    {
        if(m.vert_is_boundary(vid)) return true;
    }
    return false;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
bool vector_contains_value(std::vector<uint> &v, const uint value)
{
    bool contains=false;
    for(int i=0;i<v.size();++i)
    {
        if(v[i]==value)
        {
            contains=true;
        }
    }
    return contains;

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
        err/=2;
        break;
    }
    case 1:

    {
        norm1=v1.length();
        norm2=v2.length();
        diff_norm=abs(norm1-norm2);
        double tmp=(norm1 +norm2)/2;
        err=(tmp>treshold)? diff_norm/tmp  : treshold;
        err/=2;

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


    x=(x-md)/(Md-md);


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
DrawableVectorField  compute_ground_truth (const DrawableTrimesh<> &m, const double a, const double b, const double c, const int mode, const int method)
{
    DrawableVectorField V;
    int N=0;
    std::vector<vec3d> coords;
    double A=a*10;
    double B=b/100;
    double C=c/100;
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
        for(int i=0;i<N;++i)
        {
            vec3d pos=coords[i];
            vec3d grad=vec3d(A*B*cos(B*pos[0])*cos(B*pos[1]),-A*B*sin(B*pos[0])*sin(B*pos[1]),0);
            V.set(i,grad);
        }
        break;
    case 1:
        for (int i=0;i<N;++i)
        {
            vec3d pos=coords[i];
            vec3d grad=vec3d(2*B*(pos[0]-50),2*C*(pos[1]-50),0);
            V.set(i,grad);
        }
        break;

    case 2:
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
    double C=c/100;
    std::vector<double> data(Nv);
    switch (mode)
    {
    case 0:

        for(int i=0;i<Nv;++i)
        {
            vec3d w=v[i];
            double val=A*sin(B*w[0])*cos(B*w[1]);
            data[i]=val;
        }
        break;

    case 1:
        for(int i=0;i<Nv;++i)
        {
            vec3d w=v[i];
            double val=B*pow(w[0]-50,2)+C*pow(w[1]-50,2);
            data[i]=val;
        }
        break;

    case 2:

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

    std::vector<Eigen::ColPivHouseholderQR<Eigen::MatrixXd> > M;
    std::vector<Eigen::MatrixXd> RHS;
    std::vector<std::vector<uint>> nbrs;
    c=c.RED();
    switch(mode)
    {
    case 0:
    {
        G=gradient_matrix(m,1);
        V=DrawableVectorField(m,true);
        V=G*f;
        V.set_arrow_color(c);
    }

        break;
    case 1:
    {
        G=build_matrix_for_AGS(m);
        V=DrawableVectorField(m,false);
        V=G*f;

        V.set_arrow_color(c);
    }
        break;
    case 2:
    {
        build_matrix_for_LSDD(m,M,RHS,nbrs);

        solve_for_LSDD(m,V,M,RHS,f,nbrs);

        V.set_arrow_color(c);


    }
        break;
    case 3:
    {
        build_matrix_for_LR(m,M,RHS,nbrs);

        solve_for_LR(m,V,M,RHS,f,nbrs);

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
                avg=sqrt(avg);

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
                avg=sqrt(avg);

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
                avg=sqrt(avg);

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

        if(poly_has_vert_on_boundary(m,pid))
            m_grid.vert_data(vid).marked=true;


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




    }

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
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Eigen::SparseMatrix<double> build_matrix_for_AGS(const DrawableTrimesh<> &m)
{
    Eigen::SparseMatrix<double> G(m.num_verts()*3, m.num_verts());
    std::vector<Entry> entries;
    DrawableVectorField V=DrawableVectorField(m,false);

    uint prev;
    uint next;
    uint curr;

    vec3d c_next;
    vec3d c_prev;
    vec3d c_vert;

    for(uint vid=0; vid<m.num_verts(); ++vid)
    {

        std::vector<uint> nbr=m.vert_ordered_vert_ring(vid);
        int s=nbr.size();
        std::vector<std::pair<uint,vec3d>> vert_contr;
        double area=0.f;

        for(uint j=0;j<s;++j)
        {
            c_vert=vec3d(0,0,0);

            if(m.vert_is_boundary(vid) && j==0)
            {
                prev=vid;
                curr=nbr[j];
                next=nbr[j+1];

                uint pid=m.poly_id({prev,curr,next});
                area+=m.poly_area(pid);
                vec3d n=m.poly_data(pid).normal;

                c_prev=(m.vert(curr)-m.vert(prev)).cross(n);
                c_next=(m.vert(next)-m.vert(curr)).cross(n);

                c_vert+=c_prev;

            }else if(m.vert_is_boundary(vid) && j==s-1)
            {
                prev=nbr[j-1];
                curr=nbr[j];
                next=vid;

                uint pid=m.poly_id({prev,curr,next});
                vec3d n=m.poly_data(pid).normal;

                c_prev=(m.vert(curr)-m.vert(prev)).cross(n);
                c_next=(m.vert(next)-m.vert(curr)).cross(n);
                c_vert+=c_next;
            }
            else
            {
                prev=(j==0)?nbr[s-1]:nbr[j-1];
                curr=nbr[j];
                next=(j==s-1)? nbr[0]:nbr[j+1];

                uint pid0=m.poly_id({prev,vid,curr});
                uint pid1=m.poly_id({curr,vid,next});

                vec3d n0=m.poly_data(pid0).normal;
                vec3d n1=m.poly_data(pid1).normal;

                area+=m.poly_area(pid1);

                c_prev=(m.vert(curr)-m.vert(prev)).cross(n0);
                c_next=(m.vert(next)-m.vert(curr)).cross(n1);
            }

            vert_contr.push_back(std::make_pair(curr, (c_prev+c_next)/2));
            if(c_vert.length()>0)
                vert_contr.push_back(std::make_pair(vid, c_vert/2));
        }
        // note: Eigen::setFromTriplets will take care of summing contributs w.r.t. multiple polys

        uint row = vid * 3;
        for(auto c : vert_contr)
        {
            entries.push_back(Entry(row  , c.first, c.second.x()/area));
            entries.push_back(Entry(row+1, c.first, c.second.y()/area));
            entries.push_back(Entry(row+2, c.first, c.second.z()/area));
        }
    }

    G.setFromTriplets(entries.begin(), entries.end());

    return G;
}


//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void build_matrix_for_LSDD(DrawableTrimesh<> &m, std::vector<Eigen::ColPivHouseholderQR<MatrixXd> > &MFact, std::vector<Eigen::MatrixXd> &RHS, std::vector<std::vector<uint> > &nbrs)
{
    int Nv=m.num_verts();
    int M=0;
    double wgt;
    std::vector<uint> nbr;
    vec3d grad;
    DrawableVectorField V=DrawableVectorField(m,false);
    Eigen::VectorXd X(2);
    Eigen::VectorXd b(2);
    std::vector<double> values;
    for (int i=0;i<Nv;++i)
    {
        values.clear();
        nbr=m.vert_ordered_vert_ring(i);


        if(nbr.size()<3)
        {
            std::set<uint> two_ring=m.vert_n_ring(i,2);
            nbr.clear();
            for(uint k : two_ring)
            {

                nbr.push_back(k);

            }
        }
        int s=nbr.size();
        nbrs.push_back(nbr);
        Eigen::MatrixXd B(s,2);
        Eigen::VectorXd sigma(s);
        Eigen::DiagonalMatrix<double,Eigen::Dynamic> W(s);
        Eigen::MatrixXd Rhs;
        for(int j=0;j<s;++j)
        {
            vec3d vij=m.vert(nbr[j])-m.vert(i);
            wgt=1/vij.length_squared();
            B(j,0)=vij[0];
            B(j,1)=vij[1];
            W.diagonal()[j]=wgt;
        }

        Eigen::MatrixXd A;
        Eigen::MatrixXd Bt=Transpose<Eigen::MatrixXd>(B);

        A=Bt*W*B;
        Rhs=Bt*W;
        Eigen::ColPivHouseholderQR<MatrixXd> dec(A);
        MFact.push_back(dec);
        RHS.push_back(Rhs);

    }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void build_matrix_for_LR(DrawableTrimesh<> &m, std::vector<Eigen::ColPivHouseholderQR<Eigen::MatrixXd> > &M, std::vector<Eigen::MatrixXd> &RHS,std::vector<std::vector<uint> > &nbrs)
{
    int nv=m.num_verts();
    std::vector<uint> nbr;
    Eigen::VectorXd X(6);


    double sigma=pow(m.edge_avg_length(),2);
    double factor=1/sigma*sqrt(2*M_PI);

    for (int i=0;i<nv;++i)
    {
        nbr=m.vert_ordered_vert_ring(i);

        nbr.push_back(i);

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

        int s=nbr.size();
        nbrs.push_back(nbr);
        Eigen::MatrixXd coeff(s,6);

        Eigen::DiagonalMatrix<double,Eigen::Dynamic> W(s);

        for (int j=0;j<s;++j)
        {
            vec3d pos;
            double wgt;

            pos=m.vert(nbr[j]);
            double d=(pos-m.vert(i)).length_squared();
            d/=sigma;
            wgt=exp(-d)*factor;
            W.diagonal()[j]=wgt;



            coeff(j,0)=pow(pos[0],2);
            coeff(j,1)=pos[0]*pos[1];
            coeff(j,2)=pow(pos[1],2);
            coeff(j,3)=pos[0];
            coeff(j,4)=pos[1];
            coeff(j,5)=1;

        }
        Eigen::MatrixXd coeffT=Transpose<Eigen::MatrixXd>(coeff);
        Eigen::MatrixXd A=coeffT*W*coeff;
        Eigen::MatrixXd Rhs=coeffT*W;


        Eigen::ColPivHouseholderQR<Eigen::MatrixXd>dec(A);


        M.push_back(dec);
        RHS.push_back(Rhs);
    }
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void solve_for_LSDD(DrawableTrimesh<> &m, DrawableVectorField &V, std::vector<Eigen::ColPivHouseholderQR<MatrixXd> > &M, std::vector<MatrixXd> &RHS, ScalarField & f, std::vector<std::vector<uint> > &nbrs)
{
    V=DrawableVectorField(m,false);
    Eigen::VectorXd X(2);
    Eigen::VectorXd b(2);
    vec3d grad;


    for(int i=0;i<m.num_verts();++i)
    {
        Eigen::MatrixXd Rhs=RHS[i];

        std::vector<uint> nbr=nbrs[i];
        int p=nbr.size();

        Eigen::VectorXd tmp(p);

        for(int j=0;j<p;++j)
        {
            tmp(j)=f[nbr[j]]-f[i];
        }
        b=Rhs*tmp;


        X=M[i].solve(b);

        grad=vec3d(X(0),X(1),0);

        V.set(i,grad);


    }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void solve_for_LR(DrawableTrimesh<> &m, DrawableVectorField &V, std::vector<Eigen::ColPivHouseholderQR<Eigen::MatrixXd> > &M, std::vector<MatrixXd> RHS, ScalarField &f, std::vector<std::vector<uint> > &nbrs)
{
    V=DrawableVectorField(m,false);
    Eigen::VectorXd X(6);
    Eigen::VectorXd b(6);
    vec3d grad;
    vec3d vert;

    for(int i=0;i<RHS.size();++i)
    {
        vert=m.vert(i);
        Eigen::MatrixXd AtW=RHS[i];
        int p=nbrs[i].size();
        Eigen::VectorXd tmp(p);
        for(int j=0;j<nbrs[i].size();++j)
        {
            tmp(j)=f[nbrs[i][j]];
        }


        b=AtW*tmp;
        X=M[i].solve(b);
        grad= vec3d (2*X(0)*vert.x()+X(1)*vert.y()+X(3),2*X(2)*vert.y()+X(1)*vert.x()+X(4),0);
        V.set(i,grad);

    }


}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
