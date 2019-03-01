#include "computations.h"

#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include<Eigen/SparseQR>
#include<Eigen/OrderingMethods>
#include<iostream>
#include <cinolib/symbols.h>
#include <cinolib/color.h>
#include <cinolib/gradient.h>
#include <cinolib/point_inside_mesh.h>


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
    double angle;
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
std::vector<double> error_for_hist(const DrawableVectorField &GT, const DrawableVectorField &V, const DrawableTrimesh<> &m, const int method, const int type_of_vertices, int mode)
{
    int N=m.num_verts();
    std::vector<double> err;

    switch(type_of_vertices)
    {
    case 0:
    {

        for(int i=0;i<N;++i)
        {
            err.push_back(relative_error(GT.vec_at(i),V.vec_at(i),mode));
        }

    } break;
    case 1:
    {

        for(int i=0;i<N;++i)
        {
            if(!m.vert_is_boundary(i))
                err.push_back(relative_error(GT.vec_at(i),V.vec_at(i),mode));
        }

    } break;
    case 2:
    { for(int i=0;i<N;++i)
        {
            if(m.vert_is_boundary(i))
                err.push_back(relative_error(GT.vec_at(i),V.vec_at(i),mode));

        }

    }break;

    }
    return err;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
std::vector<double> compute_for_hist(DrawableTrimesh<> &m, DrawableTrimesh<> &grid, const int method, const double a, const double b, const int tri_type, const int N, ScalarField &F, const double anisotropy, int without_boundary, int mode, int weight)
{
    DrawableVectorField V;
    DrawableVectorField Vgrid;
    DrawableVectorField GT;


    GT=compute_ground_truth(grid,method,a,b,10,2);
    F=get_scalar_field(m,a,b,2,0,0);
    V=compute_field(m,F,method,weight);


    bring_the_field_inside(m,grid,V,Vgrid,method);
    std::vector<double> err=error_for_hist(GT,Vgrid,grid,method,without_boundary,mode);

    return err;

}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double scale_function(double x,const double scale_factor)
{

    //x=(x-md)/(Md-md)*(Mr-mr)+mr;
    double tmp=x;
    if(tmp>0.5)
    {
        x=tmp*(1+scale_factor*pow(tmp-0.5,1.5));
    }else
    {
        x=tmp*(1-scale_factor*pow(-tmp+0.5,1.5));
    }

    if(x>1)
    {
        x=1;
    }else if(x<0)
    {
        x=0;
    }



    return x;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
uint closest_vertex (const vec3d & p, const DrawableTrimesh<> *m)
{
    std::vector<std::pair<double,uint>> closest;
    for(uint vid=0; vid<m->num_verts(); ++vid) closest.push_back(std::make_pair(m->vert(vid).dist(p),vid));
    std::sort(closest.begin(), closest.end());
    return closest.front().second;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void find_max_min_values (const DrawableVectorField &V, vec3d &max, vec3d &min)
{
    int N=V.rows()/3;
    max=vec3d(0,0,0);
    min=vec3d(0,0,0);
    for(int i=0;i<N;++i)
    {
        max.x()=std::max(V.vec_at(i).x(),max.x());
        max.y()=std::max(V.vec_at(i).y(),max.y());
        max.z()=std::max(V.vec_at(i).z(),max.z());
        min.x()=std::min(V.vec_at(i).x(),min.x());
        min.y()=std::min(V.vec_at(i).y(),min.y());
        min.z()=std::min(V.vec_at(i).z(),min.z());
    }
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
DrawableVectorField from_f2v(DrawableVectorField & W, const DrawableTrimesh<> & m,int weight)

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

        provv=vec3d(0,0,0);
        avg=vec3d(0,0,0);
        vec3d pos=m.vert(i);
        double wgt=0;
        switch(weight)
        {
        case 0:
        {
            for(uint pid : m.adj_v2p(i))
            {
                provv=W.vec_at(pid)*(m.poly_area(pid));
                wgt+=m.poly_area(pid);
                avg+=provv;

            }
        }
            break;
        case 1:
        {
            for(uint pid : m.adj_v2p(i))
            {
                provv=W.vec_at(pid)*(1/(m.poly_centroid(pid)-pos).length_squared());
                wgt+=1/(m.poly_centroid(pid)-pos).length_squared();
                avg+=provv;

            }
        }
            break;
        case 2:
        {
            for(uint pid : m.adj_v2p(i))
            {
                uint eid=m.edge_opposite_to(pid,i);
                vec3d v0=m.edge_vert(eid,0)-pos;
                vec3d v1=m.edge_vert(eid,1)-pos;
                double angle=v0.angle_rad(v1);
                provv=W.vec_at(pid)*angle;
                wgt+=angle;
                avg+=provv;

            }
        }

        }

        avg/=wgt;
        Wv.set(i,avg);
    }
    return Wv;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
DrawableVectorField  compute_ground_truth (const DrawableTrimesh<> &m, const int method, const double a, const double b, const double c, const int mode)
{
    DrawableVectorField V;
    V=DrawableVectorField(m,false);
    int N=m.num_verts();
    std::vector<vec3d> coords=m.vector_verts();

    double A=a*100;
    double B=b/100;
    //     if(method==0)
    //    {
    //        V=DrawableVectorField(m,true);
    //        N=m.num_polys();
    //        for (int j=0;j<N;++j)
    //        {
    //            coords.push_back(m.poly_centroid(j));
    //        }
    //    }else
    //    {

    //    V=DrawableVectorField(m,false);
    //    N=m.num_verts();
    //    coords=m.vector_verts();
    //    }

    cinolib::Color col;
    col=col.GREEN();
    switch(mode)
    {

    case 0:
        for (int i=0;i<N;++i)
        {

            vec3d pos=coords[i];
            vec3d grad=vec3d((a*sin(3*b*pos[0])+a*3*b*pos[0]*cos(3*b*pos[0]))*sin(3*b*pos[1]*pos[1]),
                    6*b*a*pos[0]*pos[1]*sin(3*b*pos[0])*cos(3*b*pos[1]*pos[1]),0);
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
ScalarField get_scalar_field(const DrawableTrimesh<> & m, const double a, const double b, const int method,int noise,double k)
{
    std::vector<vec3d> v=m.vector_verts();
    int Nv= m.num_verts();
    double A=a*100;
    double B=b/100;

    double c=10;
    std::vector<double> data(Nv);
    switch (method)
    {
    case 0:

        for(int i=0;i<Nv;++i)
        {
            vec3d w=v[i];
            double val=a*w[0]*sin(3*b*w[0])*sin(3*b*w[1]*w[1]);
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
            double epsilon=rand() % 21 -10;//random number in the range [-10,10];
            epsilon/=100;

            double val= (noise==0) ? A*sin(B*w[0])*cos(B*w[1]) : A*(sin(B*w[0])*cos(B*w[1])+k*epsilon);

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

DrawableVectorField compute_field(DrawableTrimesh<> & m, ScalarField & f, const int method, int weight)
{

    DrawableVectorField V;
    Eigen::SparseMatrix<double> G;
    cinolib::Color c;
    c=c.BLACK();
    switch(method)
    {
    case 0:


        V=compute_PCE(m,f);
        V.set_arrow_color(c);


        break;
    case 1:

        V=compute_AGS(m,f,2);
        V.set_arrow_color(c);

        break;
    case 2:

    {

        V=compute_FEM(m,f);
        V.set_arrow_color(c);
    }


        break;
    case 3:
    {
        V=compute_quadratic_regression(m,f);
        V.set_arrow_color(c);
    }

        break;
    case 4:
    {
        V=compute_quadratic_regression_centroids(m,f);
        V.set_arrow_color(c);
    }
        break;
    }
    return V;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double estimate_MSE(const DrawableVectorField & GT, const DrawableVectorField & V, const DrawableTrimesh<> &m,const int method,const int type_of_vertices,int mode)
{

    double err;
    double max=0;
    double min=0;
    double avg=0;
    int count=0;
    int N=m.num_verts();


    //   if(method!=0)
    //    {


    switch(type_of_vertices)
    {
    case 0:
    {
        for(int i=0;i<N;++i)
        {
            err=relative_error(GT.vec_at(i),V.vec_at(i),mode);
            avg+=err;
            max=std::max(err,max);
            min=std::min(err,min);
        }
        avg/=N;
    } break;
    case 1:
    {
        for(int i=0;i<N;++i)
        {
            if(m.vert_is_boundary(i) || m.vert_data(i).marked)
            {++count;}
            else
            {
                err=relative_error(GT.vec_at(i),V.vec_at(i),mode);
                avg+=err;

                max=std::max(err,max);
                min=std::min(err,min);

            }

        }
        avg/=(N-count);
    } break;
    case 2:
    {
        for(int i=0;i<N;++i)
        {

            if(m.vert_is_boundary(i) )//|| m.vert_data(i).marked)
            {
                err=relative_error(GT.vec_at(i),V.vec_at(i),mode);
                avg+=err;
                ++count;
                max=std::max(err,max);
                min=std::min(err,min);}
            else
            {
                continue;

            }
        }
        avg=avg/count;
    }break;

    }
    if(mode==2)
    {
        avg=sqrt(avg);
    }
    std::cout<<"Err MAX="<<max<<std::endl;
    std::cout<<"Err MIN="<<min<<std::endl;
    std::cout<<"MEAN Err="<<avg<<std::endl;
    //   }
    //    else
    //    {
    //        avg=dual_error(m,GT,V,mode);
    //    }


    return avg;

}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double dual_error(const DrawableTrimesh<> &m, const DrawableVectorField &GT, const DrawableVectorField &V, const int mode)
{
    std::vector<vec3d>             dual_verts;
    std::vector<std::vector<uint>> dual_faces;
    dual_mesh(m, dual_verts,dual_faces,true);
    DrawablePolygonmesh<> dual_m;
    dual_m=DrawablePolygonmesh<>(dual_verts,dual_faces);


    double err;
    double max=0;
    double min=0;
    double avg=0;
    int count=0;
    int M=m.num_polys();

    for(int i=0;i<dual_m.num_verts();++i)
    {
        if(i<M)
        {
            err=relative_error(GT.vec_at(i),V.vec_at(i),mode);
            avg+=err;
            ++count;

            max=std::max(err,max);
            min=std::min(err,min);

        }
    }

    std::cout<<"Err MAX="<<max<<std::endl;
    std::cout<<"Err MIN="<<min<<std::endl;
    std::cout<<"MEAN Err="<<avg/count<<std::endl;
    avg=avg/count;
    if(mode==2)
    {
        avg=sqrt(avg);
    }
    return avg;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
DrawableVectorField compute_quadratic_regression(const DrawableTrimesh<> &m, const ScalarField & f)
{
    int nv=m.num_verts();
    std::vector<uint> nbr;
    vec3d vert;
    Eigen::VectorXd X(6);
    DrawableVectorField V=DrawableVectorField(m,false);
    std::vector<imaginary_vertex> nbr_aux;
    double sigma=pow(m.edge_avg_length(),2);
    double factor=1/sigma*sqrt(2*M_PI);
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

        //        if(m.vert_is_boundary(i) && handle_boundaries)
        //        {
        //            nbr_aux=nbr_for_boundaries(m,i);
        //        }

        Eigen::MatrixXd coeff(nbr.size(),6);
        Eigen::VectorXd b(nbr.size());
        Eigen::MatrixXd W(nbr.size(),nbr.size());
        int s=nbr.size();
        for (int j=0;j<nbr.size();++j)
        {
            vec3d pos;
            double wgt;

            pos=m.vert(nbr[j]);
            double d=(pos-vert).length_squared();
            d/=sigma;
            wgt=sqrt(exp(-d)*factor);
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
        X=dec.solve(B);
        vec3d sol= vec3d (2*X(0)*vert.x()+X(1)*vert.y()+X(3),2*X(2)*vert.y()+X(1)*vert.x()+X(4),0);
        V.set(i,sol);

    }

    return V;


}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
DrawableVectorField compute_quadratic_regression_centroids(DrawableTrimesh<> &m, const ScalarField & f)
{
    int nv=m.num_verts();
    std::vector<int> rank (nv);
    std::vector<uint> nbr;
    std::vector<double> values;
    std::vector<vec3d> coords=m.vector_verts();
    Eigen::VectorXd X(6);
    DrawableVectorField V=DrawableVectorField(m,false);


    for (int i=0;i<nv;++i)
    {
        nbr=m.vert_ordered_vert_ring(i);
        nbr.push_back(i);
        values.clear();
        int s=nbr.size();
        if(s<6)
        {
            for(uint j=0;j<s;++j)
            {
                for(uint pid : m.adj_v2p(nbr[j]))
                {
                    if(!m.poly_data(pid).marked && !m.poly_contains_vert(pid,i))
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

        Eigen::MatrixXd coeff(nbr.size(),6);
        Eigen::VectorXd b(nbr.size());





        for (int j=0;j<nbr.size();++j)
        {
            if(j<s)
            {
                vec3d pos=coords[nbr[j]];
                b(j)=f[nbr[j]];
                coeff(j,0)=pow(pos[0],2);
                coeff(j,1)=pos[0]*pos[1];
                coeff(j,2)=pow(pos[1],2);
                coeff(j,3)=pos[0];
                coeff(j,4)=pos[1];
                coeff(j,5)=1;
            }else
            {
                uint pid=nbr[j];
                vec3d pos=m.poly_centroid(pid);
                vec3d tmp=pos-coords[i];
                double wgt=1/tmp.length();
                b(j)=values[j-s]*wgt;
                coeff(j,0)=pow(pos[0],2)*wgt;
                coeff(j,1)=pos[0]*pos[1]*wgt;
                coeff(j,2)=pow(pos[1],2)*wgt;
                coeff(j,3)=pos[0]*wgt;
                coeff(j,4)=pos[1]*wgt;
                coeff(j,5)=1*wgt;

            }

        }
        Eigen::MatrixXd coeffT=Transpose<Eigen::MatrixXd>(coeff);
        Eigen::MatrixXd A=coeffT*coeff;
        Eigen::VectorXd B=coeffT*b;
        Eigen::ColPivHouseholderQR<MatrixXd> dec(A);
        X=dec.solve(B);
        vec3d sol= vec3d (2*X(0)*coords[i].x()+X(1)*coords[i].y()+X(3),2*X(2)*coords[i].y()+X(1)*coords[i].x()+X(4),0);
        V.set(i,sol);
    }



    return V;

}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
DrawableVectorField  compute_FEM(const DrawableTrimesh<> &m, const ScalarField & f)
{

    int Nv=m.num_verts();


    std::vector<uint> nbr;
    vec3d grad;

    DrawableVectorField V=DrawableVectorField(m,false);
    Eigen::VectorXd X(2);
    Eigen::VectorXd b(2);

    for (int i=0;i<Nv;++i)
    {

        nbr=m.vert_ordered_vert_ring(i);
        int s=nbr.size();


        Eigen::MatrixXd B(s,2);
        Eigen::VectorXd sigma(s);
        Eigen::MatrixXd W(s,s);

        for(int j=0;j<s;++j)
        {
            vec3d vij=m.vert(nbr[j])-m.vert(i);
            double wgt=1/vij.length();
            B(j,0)=vij[0]*wgt;
            B(j,1)=vij[1]*wgt;
            sigma(j)=(f[nbr[j]]-f[i])*wgt;


        }


        Eigen::Matrix2d A(2,2);

        Eigen::MatrixXd Bt=Transpose<Eigen::MatrixXd>(B);




        A=Bt*B;
        assert(A.rows()==2);
        assert(A.cols()==2);

        b=Bt*sigma;
        assert(b.rows()==2);
        assert(b.cols()==1);
        Eigen::ColPivHouseholderQR<Matrix2d> dec(A);
        X=dec.solve(b);
        grad=vec3d(X(0),X(1),0);
        V.set(i,grad);
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
ScalarField heat_map_normalization(const ScalarField &f,double scale_factor)
{

    ScalarField F=f;
    double max=0;
    double min=0;
    find_max_min_values(f,max,min);
    double delta=max-min;


    for(int i=0; i<f.rows(); ++i)
    {


        F[i]=(f[i]-min)/delta;
        F[i] =scale_function(F[i],scale_factor);

    }


    return F;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
std::vector<double> circum_radius(const DrawableTrimesh<> &m)
{
    int N=m.num_polys();
    std::vector<double> radi(N);
    for (int i=0;i<N;++i)
    {
        double L=1;
        for(uint eid : m.adj_p2e(i))
        {
            L*=m.edge_length(eid);
        }
        radi[i]=L/(4*m.poly_area(i));
    }
    return radi;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
std::vector<double> in_radius(const DrawableTrimesh<> &m)
{
    int N=m.num_polys();
    std::vector<double> radi(N);
    for (int i=0;i<N;++i)
    {
        double L=0;
        for(uint eid : m.adj_p2e(i))
        {
            L+=m.edge_length(eid)/2;
        }
        radi[i]=m.poly_area(i)/L;
    }

    return radi;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
std::vector<double> radii_ratio(const DrawableTrimesh<> &m, const std::vector<double> &R, std::vector<double> &r )
{
    int N=m.num_verts();
    std::vector<double> ratios;

    for(int i=0;i<N;++i)
    {
        if(m.vert_is_boundary(i))
        {
            continue;
        }
        else
        {
            double avg=0;
            int count=0;
            for(uint k : m.adj_v2p(i))
            {
                avg+=2*r[k]/R[k];
                ++count;
            }
            ratios.push_back(avg/count);
        }

    }
    return ratios;

}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double correlation_coefficient(const std::vector<double> &X,const std::vector<double> &Y)
{
    double sum_X=0;
    double sum_Y=0;
    double sum_XY=0;
    double squareSum_X=0;
    double squareSum_Y=0;
    int N=X.size();

    for(int i=0;i<N;++i)
    {
        sum_X+=X[i];
        sum_Y+=Y[i];
        sum_XY+=X[i]*Y[i];
        squareSum_X+=X[i]*X[i];
        squareSum_Y+=Y[i]*Y[i];


    }
    double corr=(N*sum_XY-sum_X*sum_Y)/sqrt((N*squareSum_X-pow(sum_X,2))*(N*squareSum_Y-pow(sum_Y,2)));
    return corr;

}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double average_neighborhood_area(const DrawableTrimesh<> &m)
{
    int N=m.num_polys();
    double avg=0;
    for(int i=0;i<N;++i)
    {
        avg+=m.poly_area(i);
    }
    avg=avg/N;
    return avg;
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

void bring_the_field_inside(const DrawableTrimesh<> &m,DrawableTrimesh<> &m_grid,DrawableVectorField &V,DrawableVectorField &W,const int method)
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
        double area = m.poly_area(pid)*2; // (2 is the average term : two verts for each edge)
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
DrawableVectorField compute_AGS(const DrawableTrimesh<> & m, ScalarField &f, int weight)
{
    DrawableVectorField W=DrawableVectorField(m,true);
    DrawableVectorField V=DrawableVectorField(m,false);
    W=compute_PCE(m,f);
    V=from_f2v(W,m,weight);

    return V;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
DrawableVectorField GG_on_verts(const DrawableTrimesh<> &m, ScalarField &f)
{
    Eigen::SparseMatrix<double> G(m.num_verts()*3, m.num_verts());
    std::vector<Entry> entries;
    DrawableVectorField V=DrawableVectorField(m,false);
    for(uint vid=0; vid<m.num_verts(); ++vid)
    {

        std::vector<uint> nbr=m.vert_ordered_vert_ring(vid);
        int s=nbr.size();
        std::vector<std::pair<uint,vec3d>> vert_contr;
        double area=0.f;

        for(uint j=0;j<s;++j)
        {
            //            if(!m.vert_is_boundary(i))

            //            {
            //                uint prev=(j==0)?nbr[s-1]:nbr[j-1];
            //                uint curr=nbr[j];
            //                uint next=(j==s-1)? nbr[0]:nbr[j+1];


            //                vec3d v0=m.vert(prev)-m.vert(vid);
            //                vec3d v1=m.vert(curr)-m.vert(vid);
            //                vec3d v2=m.vert(next)-m.vert(vid);

            //                vec3d n0=v0.cross(v1);
            //                n0.normalize();
            //                vec3d n1=v1.cross(v2);

            //                area+=0.5*n1.length();
            //                n1.normalize();

            //                vec3d c_prev=(m.vert(curr)-m.vert(prev)).cross(n0);
            //                vec3d c_next=(m.vert(next)-m.vert(curr)).cross(n1);

            //                vert_contr.push_back(std::make_pair(curr, (c_prev+c_next)/2));
            //            }


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
    V=G*f;
    return V;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void build_matrix_for_LSDD(DrawableTrimesh<> &m, std::vector<Eigen::ColPivHouseholderQR<Matrix2d> > &MFact, std::vector<Eigen::MatrixXd> &RHS, std::vector<std::vector<uint> > &nbrs)
{
    int Nv=m.num_verts();
    int M=0;
    double wgt;
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

                vec3d vij=coords[nbr[j]]-coords[i];
                wgt=1/vij.length_squared();
                B(j,0)=vij[0];
                B(j,1)=vij[1];

                W.diagonal()[j]=wgt;


        }

        Eigen::MatrixXd A;
        Eigen::MatrixXd Bt=Transpose<Eigen::MatrixXd>(B);

        A=Bt*W*B;
        Rhs=Bt*W;
        Eigen::ColPivHouseholderQR<Matrix2d> dec(A);
        MFact.push_back(dec);
        RHS.push_back(Rhs);



//        std::cout<<"before"<<std::endl;
//        std::cout<<dec.absDeterminant()<<std::endl;







    }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void build_matrix_for_LR(DrawableTrimesh<> &m, std::vector<Eigen::ColPivHouseholderQR<Eigen::MatrixXd> > &M, std::vector<Eigen::MatrixXd> &RHS,std::vector<std::vector<uint> > &nbrs)
{
    int nv=m.num_verts();
    std::vector<uint> nbr;
    vec3d vert;
    Eigen::VectorXd X(6);

    std::vector<imaginary_vertex> nbr_aux;
    double sigma=pow(m.edge_avg_length(),2);
    double factor=1/sigma*sqrt(2*M_PI);

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

        int s=nbr.size();
        nbrs.push_back(nbr);
        Eigen::MatrixXd coeff(s,6);

       Eigen::DiagonalMatrix<double,Eigen::Dynamic> W(s);

        for (int j=0;j<s;++j)
        {
            vec3d pos;
            double wgt;

            pos=m.vert(nbr[j]);
            double d=(pos-vert).length_squared();
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
void solve_for_LSDD(DrawableTrimesh<> &m, DrawableVectorField &V, std::vector<Eigen::ColPivHouseholderQR<Matrix2d> > &M, std::vector<MatrixXd> &RHS, ScalarField & f, std::vector<std::vector<uint> > &nbrs, std::chrono::duration<double> time_precom, std::chrono::duration<double> time_estimation)
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
void solve_for_LR(DrawableTrimesh<> &m,DrawableVectorField &V, std::vector<Eigen::ColPivHouseholderQR<Eigen::MatrixXd> > &M, std::vector<MatrixXd> RHS, ScalarField &f, std::vector<std::vector<uint> > &nbrs, std::chrono::duration<double> time_precom, std::chrono::duration<double> time_estimation)
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
