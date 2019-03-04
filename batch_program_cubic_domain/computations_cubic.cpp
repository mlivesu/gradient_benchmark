#include "computations_cubic.h"

#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <cinolib/symbols.h>
#include <cinolib/color.h>
#include <iostream>
#include <cinolib/gradient.h>
#include <cinolib/point_inside_mesh.h>


using namespace cinolib;
using namespace Eigen;
typedef Eigen::Triplet<double> Entry;
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
bool vector_contains_value(const std::vector<uint> &v, const uint value)
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
bool poly_has_vert_on_srf(const DrawableTetmesh<> &m, uint pid)
{
    for(uint vid : m.adj_p2v(pid))
    {
        if(m.vert_is_on_srf(vid)) return true;
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
double scale_function(double x, double md, double Md,const double scale_factor)
{

    double amplitude=Md-md;
    double delta=amplitude/25;
    double Mdtmp=md+scale_factor*delta;

    x=(x-md)/(Mdtmp-md);

    if (x>=1){x=0.999;}
    if (x<=0){x=1e-10;}
    return x;



    return x;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
std::vector<double> error_for_hist(const DrawableVectorField &GT, const DrawableVectorField &V, const DrawableTetmesh<> &m, const int method, const int type_of_vertices, int mode)
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
            if(!m.vert_is_on_srf(i) && !m.vert_data(i).marked)
                err.push_back(relative_error(GT.vec_at(i),V.vec_at(i),mode));
        }

    } break;
    case 2:
    { for(int i=0;i<N;++i)
        {
            if(m.vert_is_on_srf(i))
                err.push_back(relative_error(GT.vec_at(i),V.vec_at(i),mode));

        }

    }break;

    }
    return err;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
std::vector<double> compute_for_hist(DrawableTetmesh<> &m, DrawableTetmesh<> &grid, const int method, const double a, const double b, const int tri_type, const int N, ScalarField &F, const double anisotropy, int without_boundary, int mode, int weight)
{
    DrawableVectorField V;
    DrawableVectorField Vgrid;
    DrawableVectorField GT;


    GT=compute_ground_truth(grid,method,a,b,10,2);
    F=get_scalar_field(m,a,b,10,2,0,0);
    V=compute_field(m,F,method,2);
    bring_the_field_inside(m,grid,V,Vgrid,method);

    std::vector<double> err=error_for_hist(GT,Vgrid,grid,method,without_boundary,mode);

    return err;

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
double solid_angle(const DrawableTetmesh<> &m, uint &pid, uint &vid)
{
    double angle=0;
    std::vector<uint> eids;
    for(uint eid : m.adj_v2e(vid))
    {
        if(m.poly_contains_edge(pid,eid))
            eids.push_back(eid);
    }

    for(uint i=0;i<eids.size();++i)
    {
        std::vector<uint> fids=m.poly_e2f(pid,eids[i]);
        angle+=m.poly_dihedral_angle(pid,fids[0],fids[1]);
    }

    angle-=M_PI;
    return angle;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
DrawableVectorField from_p2v(DrawableVectorField &W, const DrawableTetmesh<> &m, int weight)

{
    int Nv=m.num_verts();
    std::vector<vec3d> grad_v(Nv);
    vec3d provv=vec3d(0,0,0);
    vec3d avg=vec3d(0,0,0);
    double wgt;
    DrawableVectorField  Wv;
    Wv=DrawableVectorField (m,false);


    for (uint vid=0;vid<Nv;++vid)
    {
        provv=vec3d(0,0,0);
        avg=vec3d(0,0,0);
        wgt=0;
        vec3d pos=m.vert(vid);

        switch(weight)
        {
        case 0:
        {
            for(uint pid : m.adj_v2p(vid))
            {
                provv=W.vec_at(pid)*(m.poly_volume(pid));
                wgt+=m.poly_volume(pid);
                avg+=provv;

            }
        }
            break;
        case 1:
        {
            for(uint pid : m.adj_v2p(vid))
            {
                provv=W.vec_at(pid)*(1/(m.poly_centroid(pid)-pos).length_squared());
                wgt+=1/(m.poly_centroid(pid)-pos).length_squared();
                avg+=provv;

            }
            break;
        }
        case 2:
        {
            for(uint pid : m.adj_v2p(vid))
            {
                double angle=solid_angle(m,pid,vid);
                provv=W.vec_at(pid)*angle;
                wgt+=angle;
                avg+=provv;

            }
            break;
        }
        }

        avg/=wgt;
        Wv.set(vid,avg);
    }
    return Wv;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
DrawableVectorField  compute_ground_truth (const DrawableTetmesh<> &m, const int method, const double a,const double b, const double c, const int mode)
{
    DrawableVectorField V;
    int N=0;
    std::vector<vec3d> coords;
    //    if(method==0)
    //    {
    //        V=DrawableVectorField(m,true);
    //        N=m.num_polys();
    //        for (int j=0;j<N;++j)
    //        {
    //            coords.push_back(m.poly_centroid(j));
    //        }
    //    }else
    //    {

    //        V=DrawableVectorField(m,false);
    //        N=m.num_verts();
    //        coords=m.vector_verts();
    //    }


    V=DrawableVectorField(m,false);
    N=m.num_verts();
    coords=m.vector_verts();
    cinolib::Color col;
    col=col.GREEN();
    switch(mode)
    {

    case 0:
        for (int i=0;i<N;++i)
        {
            double B=b/10;
            double A=a/10;
            vec3d pos=coords[i];
            double tmp=A*sin(2*M_PI*B*cos(M_PI*sqrt(pow(pos[0],2)+pow(pos[1],2))/2))*pow(M_PI,2)*B*sin(M_PI*sqrt(pow(pos[0],2)+pow(pos[1],2))/2)/((sqrt(pow(pos[0],2)+pow(pos[1],2)))*(2+2*A));
            vec3d grad=vec3d(tmp*pos[0],tmp*pos[1],-cos(M_PI*pos[2]/2)*M_PI/(4*(1+A)));
            V.set(i,grad);


        }
        break;
    case 1:
        for (int i=0;i<N;++i)
        {

            vec3d pos=coords[i];
            vec3d grad=vec3d(2*b*(pos[0]-0.5),2*c*(pos[1]-0.5),2*b*(pos[2]-0.5));
            V.set(i,grad);


        }
        break;
    case 2:
        for(int i=0;i<N;++i)
        {

            double A=a*100;
            double B=b/100;
            vec3d pos=coords[i];
            vec3d grad=vec3d(A*B*cos(B*pos[0])*cos(B*pos[1])*sin(B*pos[2]),-A*B*sin(B*pos[0])*sin(B*pos[1])*sin(B*pos[2]),A*B*sin(B*pos[0])*cos(B*pos[1])*cos(B*pos[2]));
            V.set(i,grad);

        }
        break;
    case 3:
        for (int i=0;i<N;++i)
        {
            vec3d grad=vec3d(1,1,1);
            V.set(i,grad);




        }
        break;
    }


    V.set_arrow_color(col);
    return V;

}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
ScalarField get_scalar_field(const DrawableTetmesh<> &m, const double a,const double b,const double c, const int method,int noise,double k)
{

    std::vector<vec3d> v=m.vector_verts();
    int Nv= m.num_verts();



    std::vector<double> data(Nv);
    switch (method)
    {
    case 0:

        for(int i=0;i<Nv;++i)
        {
            double B=b/10;
            double A=a/10;
            vec3d w=v[i];
            double val=(1-sin(M_PI*w[2]/2)+A*(1+cos(2*M_PI*B*cos(M_PI*sqrt(pow(w[0],2)+pow(w[1],2))/2))))/(2*(1+A));
            data[i]=val;


        }
        break;
    case 1:
        for(int i=0;i<Nv;++i)
        {

            vec3d w=v[i];
            double val=b*pow(w[0]-0.5,2)+c*pow(w[1]-0.5,2)+b*pow(w[2]-0.5,2);
            data[i]=val;


        }
        break;
    case 2:
        for(int i=0;i<Nv;++i)
        {
            vec3d w=v[i];
            double A=a*100;
            double B=b/100;
            double epsilon=rand() % 21 -10;
            epsilon/=100;

            double val= (noise==0) ? A*sin(B*w[0])*cos(B*w[1])*sin(B*w[2]) : A*(sin(B*w[0])*cos(B*w[1])*sin(B*w[2])+k*epsilon);

            data[i]=val;
        }
        break;
    case 3:

        for(int i=0;i<Nv;++i)
        {

            vec3d w=v[i];
            double val=w[0]+w[1]+w[2];
            data[i]=val;



        }

        break;
    }


    return ScalarField(data);
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

DrawableVectorField compute_field(const DrawableTetmesh<> & m, ScalarField & f, const int method, int weight)
{

    DrawableVectorField V;
    DrawableVectorField W;
    Eigen::SparseMatrix<double> G;
    cinolib::Color c;
    c=c.BLACK();
    switch(method)
    {
    case 0:


        G=gradient_matrix(m);
        V=DrawableVectorField(m,true);
        V=G*f;
        V.set_arrow_color(c);

        break;
    case 1:

        G=build_matrix_for_AGS(m);
        V=DrawableVectorField(m,false);
        V=G*f;
        V.set_arrow_color(c);

        break;
    case 2:

        V=compute_FEM(m,f);
        V.set_arrow_color(c);


        break;
    case 3:

        V=compute_quadratic_regression(m,f);
        V.set_arrow_color(c);

        break;
    case 4:
        V=DrawableVectorField(m,false);
        V=compute_FEM_with_centroids(m,f);
        V.set_arrow_color(c);

        break;
    case 5:
        V=DrawableVectorField(m,false);
        V=compute_quadratic_regression_with_centroids(m,f);
        V.set_arrow_color(c);

        break;
    }
    return V;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double estimate_MSE(const DrawableVectorField & GT, const DrawableVectorField & V, const DrawableTetmesh<> &grid, const int method, const int type_of_vertices, int mode)
{

    double err;
    double max=0;
    double min=1e+5;
    double avg=0;
    int count=0;

    int N=grid.num_verts();

    //    if(method!=0)
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
            if(grid.vert_is_on_srf(i) || grid.vert_data(i).marked)
            {
                ++count;
            }
            else
            {
                err=relative_error(GT.vec_at(i),V.vec_at(i),mode);
                avg+=err;
                max=std::max(err,max);
                min=std::min(err,min);
            }

        }
        avg=avg/(N-count);
    } break;
    case 2:
    { for(int i=0;i<N;++i)
        {
            if(grid.vert_is_on_srf(i))
            {
                err=relative_error(GT.vec_at(i),V.vec_at(i),mode);
                avg+=err;
                max=std::max(err,max);
                min=std::min(err,min);
            }
            else
            {
                ++count;
            }

        }
        avg=avg/(N-count);
    }break;

    }
    std::cout<<"Err MAX="<<max<<std::endl;
    std::cout<<"Err MIN="<<min<<std::endl;
    std::cout<<"MEAN Err="<<avg<<std::endl;


    if(mode==2)
    {
        avg=sqrt(avg);
    }




    //    }else
    //    {
    //        avg=dual_error(grid,GT,V,mode);
    //    }


    return avg;

}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double dual_error(const DrawableTetmesh<> &m, const DrawableVectorField &GT, const DrawableVectorField &V, const int mode)
{
    std::vector<vec3d>             dual_verts;
    std::vector<std::vector<uint>> dual_faces;
    DrawablePolyhedralmesh<> dual_m;
    dual_mesh(m, dual_m,true);

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

        }else
        {
            continue;
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
DrawableVectorField compute_quadratic_regression(const DrawableTetmesh<> &m, const ScalarField & f)
{
    int nv=m.num_verts();
    std::vector<int> rank (nv);
    std::vector<uint> nbr;
    std::vector<double> values;
    std::vector<std::pair<double,vec3d>> nbr_aux;
    std::vector<vec3d> coords=m.vector_verts();
    Eigen::VectorXd X(10);
    double sigma=pow(m.edge_avg_length(),2);
    double factor=1/sigma*sqrt(2*M_PI);
    DrawableVectorField V=DrawableVectorField(m,false);

    double count=0;

    for (int i=0;i<nv;++i)
    {
        nbr.resize(0);
        nbr_aux.resize(0);
        values.resize(0);
        vec3d vert=m.vert(i);
        for(uint vid : m.adj_v2v(i))
        {
            nbr.push_back(vid);

        }
        nbr.push_back(i);
        if(nbr.size()<9)
               {

                   for(int s=0;s<nbr.size()-1;++s)
                   {
                       if(nbr.size()>=9)
                       {
                           break;
                       }
                       else
                       {
                           for(uint vid : m.adj_v2v(nbr[s]))
                           {

                               if(vector_contains_value(nbr,vid))
                               {
                                   continue;
                               }
                               else
                               {
                                   nbr.push_back(vid);
                               }

                           }

                       }

                   }
               }
        /*if(nbr.size()<9)
        {
            std::vector<uint> picked_ones;
            for(int s=0;s<nbr.size()-1;++s)
            {
                if(nbr.size()+nbr_aux.size()>=16)
                {
                    break;
                }else
                {
                    for(uint pid : m.adj_v2p(nbr[s]))
                    {
                        if(m.poly_contains_vert(pid,i)||vector_contains_value(picked_ones,pid))
                        {
                            continue;
                        }else
                        {
                            picked_ones.push_back(pid);

                            std::vector<uint> pid_vertices=m.poly_verts_id(pid);
                            double value=0;
                            for(uint j=0;j<pid_vertices.size();++j)
                            {
                                value+=f[pid_vertices[j]];
                            }
                            nbr_aux.push_back(std::make_pair(value/4, m.poly_centroid(pid)));

                        }
                    }
                }
            }
        }*/

        int size=nbr.size()+nbr_aux.size();
        Eigen::MatrixXd coeff(size,10);
        Eigen::VectorXd b(size);

        for (int j=0;j<size;++j)
        {
            vec3d pos;
            double wgt;
            if (j<nbr.size())
            {
                pos=coords[nbr[j]];
                double d=(pos-vert).length_squared();
                d/=sigma;
                wgt=sqrt(exp(-d)*factor);

                b(j)=f[nbr[j]]*wgt;
            }else
            {
                pos=nbr_aux[j-nbr.size()].second;
                b(j)=nbr_aux[j-nbr.size()].first;
            }


            coeff(j,0)=pow(pos[0],2)*wgt;
            coeff(j,1)=pow(pos[1],2)*wgt;
            coeff(j,2)=pow(pos[2],2)*wgt;
            coeff(j,3)=pos[0]*pos[1]*wgt;
            coeff(j,4)=pos[0]*pos[2]*wgt;
            coeff(j,5)=pos[2]*pos[1]*wgt;
            coeff(j,6)=pos[0]*wgt;
            coeff(j,7)=pos[1]*wgt;
            coeff(j,8)=pos[2]*wgt;
            coeff(j,9)=1*wgt;
        }


        Eigen::MatrixXd coeffT=Transpose<Eigen::MatrixXd>(coeff);
        Eigen::MatrixXd A=coeffT*coeff;
        Eigen::VectorXd B=coeffT*b;
        Eigen::ColPivHouseholderQR<MatrixXd> dec(A);
        X=dec.solve(B);

        X.resize(10);
        vec3d sol= vec3d (2*X(0)*coords[i].x()+X(3)*coords[i].y()+X(4)*coords[i].z()+X(6),2*X(1)*coords[i].y()+X(3)*coords[i].x()+X(5)*coords[i].z()+X(7),2*X(2)*coords[i].z()+X(4)*coords[i].x()+X(5)*coords[i].y()+X(8));
        V.set(i,sol);


    }

    return V;


}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
DrawableVectorField compute_quadratic_regression_with_centroids(const DrawableTetmesh<> &m, const ScalarField & f)
{
    int nv=m.num_verts();
    std::vector<int> rank (nv);
    std::vector<std::pair<double,vec3d>> nbr;
    std::vector<double> values;
    std::vector<vec3d> coords=m.vector_verts();
    Eigen::VectorXd X(10);
    DrawableVectorField V=DrawableVectorField(m,false);




    for (int i=0;i<nv;++i)
    {
        nbr.resize(0);
        values.resize(0);
        double max=0;
        for(uint pid : m.adj_v2p(i))
        {


            std::vector<uint> pid_vertices=m.poly_verts_id(pid);
            double value=0;
            for(uint j=0;j<pid_vertices.size();++j)
            {
                value+=f[pid_vertices[j]];
            }
            nbr.push_back(std::make_pair(value/4, m.poly_centroid(pid)));



        }
        nbr.push_back(std::make_pair(f[i], coords[i]));


        if(nbr.size()<9)
        {
            std::vector<uint> picked_ones;
            for (uint vid : m.adj_v2v(i))
            {
                if(nbr.size()>=9)
                {
                    break;
                }else
                {
                    for(uint pid : m.adj_v2p(vid))
                    {
                        if(m.poly_contains_vert(pid,i))
                        {
                            continue;
                        }else
                        {
                            picked_ones.push_back(pid);

                            std::vector<uint> pid_vertices=m.poly_verts_id(pid);
                            double value=0;
                            for(uint j=0;j<pid_vertices.size();++j)
                            {
                                value+=f[pid_vertices[j]];
                            }
                            nbr.push_back(std::make_pair(value/4, m.poly_centroid(pid)));

                        }
                    }
                }



            }
        }

        Eigen::MatrixXd coeff(nbr.size(),10);
        Eigen::VectorXd b(nbr.size());


        for (int j=0;j<nbr.size();++j)
        {


            //double wgt=m.poly_volume(nbr[j].first);

            //max=std::max(max,wgt);
            vec3d pos=nbr[j].second;
            vec3d tmp=pos.operator -(coords[i]);

            double wgt=1;//1+tmp.length_squared();
            b(j)=nbr[j].first/wgt;
            coeff(j,0)=pow(pos[0],2)/wgt;
            coeff(j,1)=pow(pos[1],2)/wgt;
            coeff(j,2)=pow(pos[2],2)/wgt;
            coeff(j,3)=pos[0]*pos[1]/wgt;
            coeff(j,4)=pos[0]*pos[2]/wgt;
            coeff(j,5)=pos[2]*pos[1]/wgt;
            coeff(j,6)=pos[0]/wgt;
            coeff(j,7)=pos[1]/wgt;
            coeff(j,8)=pos[2]/wgt;
            coeff(j,9)=1/wgt;
        }


        Eigen::MatrixXd coeffT=Transpose<Eigen::MatrixXd>(coeff);
        Eigen::MatrixXd A=coeffT*coeff;
        Eigen::VectorXd B=coeffT*b;
        Eigen::ColPivHouseholderQR<MatrixXd> dec(A);
        X=dec.solve(B);
        X.resize(10);
        vec3d sol= vec3d (2*X(0)*coords[i].x()+X(3)*coords[i].y()+X(4)*coords[i].z()+X(6),2*X(1)*coords[i].y()+X(3)*coords[i].x()+X(5)*coords[i].z()+X(7),2*X(2)*coords[i].z()+X(4)*coords[i].x()+X(5)*coords[i].y()+X(8));
        V.set(i,sol);
    }




    return V;

}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
DrawableVectorField compute_FEM_with_centroids(const DrawableTetmesh<> &m, const ScalarField & f)
{


    int Nv=m.num_verts();
    int M=0;
    double delta=0;
    std::vector<vec3d> nbr;
    std::vector<double> values;

    vec3d grad;
    std::vector<vec3d> coords=m.vector_verts();
    DrawableVectorField V=DrawableVectorField(m,false);
    Eigen::Vector3d X;
    Eigen::Vector3d b;


    for (int i=0;i<Nv;++i)
    {
        nbr.resize(0);
        values.resize(0);

        for(uint pid : m.adj_v2p(i))
        {
            nbr.push_back(m.poly_centroid(pid));
            std::vector<uint> pid_vertices=m.poly_verts_id(pid);
            double value=0;
            for(uint j=0;j<pid_vertices.size();++j)
            {
                value+=f[pid_vertices[j]];
            }
            values.push_back(value/4);


        }



        M=nbr.size();
        Eigen::MatrixXd B(M,3);
        Eigen::VectorXd sigma(M);
        for(int j=0;j<M;++j)
        {
            vec3d vij=nbr[j].operator -(coords[i]);
            B(j,0)=vij[0];
            B(j,1)=vij[1];
            B(j,2)=vij[2];
            sigma(j)=(values[j]-f[i]);

            if(j==M-1)
            {
                Eigen::Matrix3d A(3,3);
                Eigen::MatrixXd Bt=Transpose<Eigen::MatrixXd>(B);
                A=Bt*B;
                Eigen::ColPivHouseholderQR<Matrix3d> dec(A);
                b=Bt*sigma;

                X=dec.solve(b);
                // X.resize(3);

                grad=vec3d(X(0),X(1),X(2));
                V.set(i,grad);


            }
        }

    }


    return V;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
DrawableVectorField compute_FEM(const DrawableTetmesh<> &m, const ScalarField & f)
{


    int Nv=m.num_verts();
    int M=0;
    double delta=0;
    std::vector<uint> nbr;

    vec3d grad;
    std::vector<vec3d> coords=m.vector_verts();
    DrawableVectorField V=DrawableVectorField(m,false);
    Eigen::Vector3d X;
    Eigen::Vector3d b;


    for (int i=0;i<Nv;++i)
    {
        nbr.resize(0);
        for(uint vid : m.adj_v2v(i))
        {
            nbr.push_back(vid);

        }



        M=nbr.size();
        Eigen::MatrixXd B(M,3);
        Eigen::VectorXd sigma(M);
        for(int j=0;j<M;++j)
        {
            vec3d vij=coords[nbr[j]].operator -(coords[i]);
            double wgt=1/vij.length();
            B(j,0)=vij[0]*wgt;
            B(j,1)=vij[1]*wgt;
            B(j,2)=vij[2]*wgt;
            sigma(j)=(f[nbr[j]]-f[i])*wgt;
        }
        Eigen::Matrix3d A(3,3);
        Eigen::MatrixXd Bt=Transpose<Eigen::MatrixXd>(B);
        A=Bt*B;
        Eigen::ColPivHouseholderQR<Matrix3d> dec(A);
        b=Bt*sigma;

        X=dec.solve(b);

        grad=vec3d(X(0),X(1),X(2));
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
ScalarField heat_map_normalization(const ScalarField &f,double min,double max,double scale_factor)
{

    ScalarField F=f;
    for(int i=0; i<f.rows(); ++i)
    {
        F[i] =scale_function(F[i],min,max,scale_factor);

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
std::vector<double> barycentric_coordinates(const vec3d &A, const vec3d &B, const vec3d &C, const vec3d &D, const vec3d &P)
{
    std::vector<double> wgts(4, 0.0);

    //Eigen::Matrix<double,4,4> M;
    //M(0,0) = A[0]; M(0,1) = A[1]; M(0,2) = A[2]; M(0,3) = 1;
    //M(1,0) = B[0]; M(1,1) = B[1]; M(1,2) = B[2]; M(1,3) = 1;
    //M(2,0) = C[0]; M(2,1) = C[1]; M(2,2) = C[2]; M(2,3) = 1;
    //M(3,0) = D[0]; M(3,1) = D[1]; M(3,2) = D[2]; M(3,3) = 1;

    Eigen::Matrix<double,4,4> M0;
    M0(0,0) = P[0]; M0(0,1) = P[1]; M0(0,2) = P[2]; M0(0,3) = 1;
    M0(1,0) = B[0]; M0(1,1) = B[1]; M0(1,2) = B[2]; M0(1,3) = 1;
    M0(2,0) = C[0]; M0(2,1) = C[1]; M0(2,2) = C[2]; M0(2,3) = 1;
    M0(3,0) = D[0]; M0(3,1) = D[1]; M0(3,2) = D[2]; M0(3,3) = 1;

    Eigen::Matrix<double,4,4> M1;
    M1(0,0) = A[0]; M1(0,1) = A[1]; M1(0,2) = A[2]; M1(0,3) = 1;
    M1(1,0) = P[0]; M1(1,1) = P[1]; M1(1,2) = P[2]; M1(1,3) = 1;
    M1(2,0) = C[0]; M1(2,1) = C[1]; M1(2,2) = C[2]; M1(2,3) = 1;
    M1(3,0) = D[0]; M1(3,1) = D[1]; M1(3,2) = D[2]; M1(3,3) = 1;

    Eigen::Matrix<double,4,4> M2;
    M2(0,0) = A[0]; M2(0,1) = A[1]; M2(0,2) = A[2]; M2(0,3) = 1;
    M2(1,0) = B[0]; M2(1,1) = B[1]; M2(1,2) = B[2]; M2(1,3) = 1;
    M2(2,0) = P[0]; M2(2,1) = P[1]; M2(2,2) = P[2]; M2(2,3) = 1;
    M2(3,0) = D[0]; M2(3,1) = D[1]; M2(3,2) = D[2]; M2(3,3) = 1;

    Eigen::Matrix<double,4,4> M3;
    M3(0,0) = A[0]; M3(0,1) = A[1]; M3(0,2) = A[2]; M3(0,3) = 1;
    M3(1,0) = B[0]; M3(1,1) = B[1]; M3(1,2) = B[2]; M3(1,3) = 1;
    M3(2,0) = C[0]; M3(2,1) = C[1]; M3(2,2) = C[2]; M3(2,3) = 1;
    M3(3,0) = P[0]; M3(3,1) = P[1]; M3(3,2) = P[2]; M3(3,3) = 1;

    //double det_M  = M.determinant();
    double det_M0 = M0.determinant();
    double det_M1 = M1.determinant();
    double det_M2 = M2.determinant();
    double det_M3 = M3.determinant();
    double sum    = det_M0 + det_M1 + det_M2 + det_M3;



    wgts[0] = det_M0/sum; assert(!std::isnan(wgts[0]));
    wgts[1] = det_M1/sum; assert(!std::isnan(wgts[1]));
    wgts[2] = det_M2/sum; assert(!std::isnan(wgts[2]));
    wgts[3] = det_M3/sum; assert(!std::isnan(wgts[3]));

    return wgts;


}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void bring_the_field_inside(const DrawableTetmesh<> &m,DrawableTetmesh<> &m_grid,DrawableVectorField &V,DrawableVectorField &W,const int method)
{
    int Nv=m_grid.num_verts();
    W=DrawableVectorField(m_grid,false);
    PointInsideMeshCache<Tetmesh<>> cache(m);

    for(uint vid=0;vid<Nv;++vid)
    {

        uint pid;
        vec3d pos=m_grid.vert(vid);
        std::vector<double> wgts;
        cache.locate(m_grid.vert(vid), pid, wgts);

        if(poly_has_vert_on_srf(m,pid))
        {
            m_grid.vert_data(vid).marked=true;
            continue;
        }else{

            if(method==0)
            {
                W.set(vid,V.vec_at(pid));
            }else
            {


                std::vector<uint>poly_vertices_id=m.poly_verts_id(pid,true);
                std::vector<vec3d> vertices_coords(poly_vertices_id.size());
                std::vector<double> bary_coords;
                for (uint i=0;i<poly_vertices_id.size();++i)
                {
                    vertices_coords[i]= m.vert(poly_vertices_id[i]);
                }

                bary_coords=barycentric_coordinates(vertices_coords[0],vertices_coords[1],vertices_coords[2],vertices_coords[3],pos);

                vec3d interpolated_value=bary_coords[0]*V.vec_at(poly_vertices_id[0])+bary_coords[1]*V.vec_at(poly_vertices_id[1])+bary_coords[2]*V.vec_at(poly_vertices_id[2])+bary_coords[3]*V.vec_at(poly_vertices_id[3]);
                W.set(vid,interpolated_value);
            }

        }


    }



}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


DrawableVectorField compute_PCE(const DrawableTetmesh<> & m, const ScalarField &f)
{
    DrawableVectorField W=DrawableVectorField(m,true);
    int Nf=m.num_faces();
    std::vector<vec3d> normals(Nf);



    for (uint fid=0;fid<Nf;++fid)
    {
        normals[fid]=m.face_data(fid).normal;
    }


    for(int pid=0; pid<m.num_polys(); ++pid)
    {
        double vol = std::max(m.poly_volume(pid), 1e-5);
        vec3d contribute(0,0,0);

        for(uint fid : m.adj_p2f(pid))
        {
            vec3d n(0,0,0);


            if (m.poly_face_is_CCW(pid,fid))
            {
                n=normals.at(fid);
            }else
            {
                n=-normals.at(fid);
            }

            double a=m.face_area(fid);
            std::vector<uint> verts=m.face_verts_id(fid);
            double value=(f[verts[0]]+f[verts[1]]+f[verts[2]])/3;
            contribute+=n*a*value;
        }
        contribute/=vol;


        W.set(pid,contribute);
    }
    return W;

}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
DrawableVectorField compute_AGS(const DrawableTetmesh<> & m, ScalarField &f, int weight)
{
    DrawableVectorField W=DrawableVectorField(m,true);
    DrawableVectorField V=DrawableVectorField(m,false);
    W=compute_PCE(m,f);
    V=from_p2v(W,m,weight);

    return V;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Eigen::SparseMatrix<double> build_matrix_for_AGS(const DrawableTetmesh<> &m)
{
    Eigen::SparseMatrix<double> G(m.num_verts()*3, m.num_verts());
    std::vector<Entry> entries;
    DrawableVectorField V=DrawableVectorField(m,false);
    for(uint vid=0; vid<m.num_verts(); ++vid)
    {
        double vol=0;
        std::vector<std::pair<uint,vec3d>> face_contr;
        double area=0.f;
        for(uint pid : m.adj_v2p(vid))
        {
            vol+=m.poly_volume(pid);
            uint f=m.poly_face_opposite_to(pid,vid);
            vec3d n=m.poly_face_normal(pid,f);
            double a=m.face_area(f);
            vec3d contribute=(n*a)/3;

            face_contr.push_back(std::make_pair(f, contribute));

            if(m.vert_is_on_srf(vid))
            {
               std::vector<uint> faces=m.poly_f2f(pid,f);
               for(uint k=0;k<faces.size();++k)
               {
                   if(m.face_is_on_srf(faces[k]))
                   {
                       n=m.poly_face_normal(pid,faces[k]);
                       a=m.face_area(faces[k]);
                       contribute=(n*a)/3;
                       face_contr.push_back(std::make_pair(faces[k], contribute));

                   }
               }

            }



        }
        uint row = vid * 3;
        for(uint s=0;s<face_contr.size();++s)
        {
             std::vector<uint> vids=m.face_verts_id(face_contr[s].first);
             for(uint j=0;j<3;++j)
             {
                 entries.push_back(Entry(row  , vids[j],face_contr[s].second.x()/vol));
                 entries.push_back(Entry(row+1, vids[j], face_contr[s].second.y()/vol));
                 entries.push_back(Entry(row+2, vids[j], face_contr[s].second.z()/vol));

             }
        }
        // note: Eigen::setFromTriplets will take care of summing contributs w.r.t. multiple polys


    }


    G.setFromTriplets(entries.begin(), entries.end());

    return G;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void build_matrix_for_LSDD(DrawableTetmesh<> &m, std::vector<Eigen::ColPivHouseholderQR<MatrixXd> > &MFact, std::vector<Eigen::MatrixXd> &RHS, std::vector<std::vector<uint> > &nbrs)
{
    int Nv=m.num_verts();
    int M=0;
    double delta=0;
    std::vector<uint> nbr;

    vec3d grad;
    std::vector<vec3d> coords=m.vector_verts();
    DrawableVectorField V=DrawableVectorField(m,false);
    Eigen::Vector3d X;
    Eigen::Vector3d b;


    for (int i=0;i<Nv;++i)
    {
        nbr.resize(0);
        for(uint vid : m.adj_v2v(i))
        {
            nbr.push_back(vid);

        }



        M=nbr.size();
        nbrs.push_back(nbr);
        Eigen::MatrixXd B(M,3);
        Eigen::VectorXd sigma(M);
        Eigen::DiagonalMatrix<double,Eigen::Dynamic> W(M);
        for(int j=0;j<M;++j)
        {
            vec3d vij=coords[nbr[j]].operator -(coords[i]);
            double wgt=1/vij.length_squared();
            W.diagonal()[j]=wgt;
            B(j,0)=vij[0];
            B(j,1)=vij[1];
            B(j,2)=vij[2];

        }
        Eigen::MatrixXd A;
        Eigen::MatrixXd Bt=Transpose<Eigen::MatrixXd>(B);
        A=Bt*W*B;
        Eigen::MatrixXd Rhs=Bt*W;
        Eigen::ColPivHouseholderQR<MatrixXd> dec(A);
        MFact.push_back(dec);
        RHS.push_back(Rhs);
    }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void build_matrix_for_LR(DrawableTetmesh<> &m, std::vector<Eigen::ColPivHouseholderQR<Eigen::MatrixXd> > &M, std::vector<Eigen::MatrixXd> &RHS,std::vector<std::vector<uint> > &nbrs)
{
    int nv=m.num_verts();
    std::vector<int> rank (nv);
    std::vector<uint> nbr;
    std::vector<double> values;
    std::vector<std::pair<double,vec3d>> nbr_aux;
    std::vector<vec3d> coords=m.vector_verts();
    Eigen::VectorXd X(10);
    double sigma=pow(m.edge_avg_length(),2);
    double factor=1/sigma*sqrt(2*M_PI);
    DrawableVectorField V=DrawableVectorField(m,false);

    double count=0;

    for (int i=0;i<nv;++i)
    {
        nbr.resize(0);
        nbr_aux.resize(0);
        values.resize(0);
        vec3d vert=m.vert(i);
        for(uint vid : m.adj_v2v(i))
        {
            nbr.push_back(vid);

        }
        nbr.push_back(i);
        if(nbr.size()<9)
               {

                   for(int s=0;s<nbr.size()-1;++s)
                   {
                       if(nbr.size()>=9)
                       {
                           break;
                       }
                       else
                       {
                           for(uint vid : m.adj_v2v(nbr[s]))
                           {

                               if(vector_contains_value(nbr,vid))
                               {
                                   continue;
                               }
                               else
                               {
                                   nbr.push_back(vid);
                               }

                           }

                       }

                   }
               }


        int size=nbr.size();
        nbrs.push_back(nbr);
        Eigen::MatrixXd coeff(size,10);
        Eigen::VectorXd b(size);
        Eigen::DiagonalMatrix<double,Eigen::Dynamic> W(size);

        for (int j=0;j<size;++j)
        {
            vec3d pos;
            double wgt;
            pos=coords[nbr[j]];
            double d=(pos-vert).length_squared();
            d/=sigma;
            wgt=exp(-d)*factor;
            W.diagonal()[j]=wgt;



            coeff(j,0)=pow(pos[0],2);
            coeff(j,1)=pow(pos[1],2);
            coeff(j,2)=pow(pos[2],2);
            coeff(j,3)=pos[0]*pos[1];
            coeff(j,4)=pos[0]*pos[2];
            coeff(j,5)=pos[2]*pos[1];
            coeff(j,6)=pos[0];
            coeff(j,7)=pos[1];
            coeff(j,8)=pos[2];
            coeff(j,9)=1;
        }


        Eigen::MatrixXd coeffT=Transpose<Eigen::MatrixXd>(coeff);
        Eigen::MatrixXd A=coeffT*W*coeff;
        Eigen::MatrixXd Rhs=coeffT*W;
        Eigen::ColPivHouseholderQR<MatrixXd> dec(A);
        M.push_back(dec);
        RHS.push_back(Rhs);
      }
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void solve_for_LSDD(DrawableTetmesh<> &m, DrawableVectorField &V, std::vector<Eigen::ColPivHouseholderQR<MatrixXd> > &M, std::vector<MatrixXd> &RHS, ScalarField & f, std::vector<std::vector<uint> > &nbrs)
{
    V=DrawableVectorField(m,false);
    Eigen::VectorXd X(3);
    Eigen::VectorXd b(3);
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

       grad=vec3d(X(0),X(1),X(2));

       V.set(i,grad);


   }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void solve_for_LR(DrawableTetmesh<> &m, DrawableVectorField &V, std::vector<Eigen::ColPivHouseholderQR<Eigen::MatrixXd> > &M, std::vector<MatrixXd> &RHS, ScalarField &f, std::vector<std::vector<uint> > &nbrs)
{
   V=DrawableVectorField(m,false);
    Eigen::VectorXd X(10);
     Eigen::VectorXd b(10);
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
      grad= vec3d (2*X(0)*m.vert(i).x()+X(3)*m.vert(i).y()+X(4)*m.vert(i).z()+X(6),2*X(1)*m.vert(i).y()+X(3)*m.vert(i).x()+X(5)*m.vert(i).z()+X(7),2*X(2)*m.vert(i).z()+X(4)*m.vert(i).x()+X(5)*m.vert(i).y()+X(8));
      V.set(i,grad);

    }



}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
