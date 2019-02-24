#include <QApplication>
#include<iostream>
#include<iomanip>
#include<computations.h>
#include<triangulations.h>
#include<functions.h>
#include <cinolib/scalar_field.h>
#include <cinolib/drawable_vector_field.h>
#include<cinolib/meshes/abstract_mesh.h>
#include <ctime>
#include <ratio>
#include <chrono>


using namespace cinolib;
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double compute(DrawableTrimesh<> &m, DrawableTrimesh<> &grid,const int method,const double a, const double b,const int tri_type, const int N, ScalarField &F,const double anisotropy,int without_boundary,int mode,int weight)
{


   // make_triangulation(m,tri_type,N,anisotropy);

    DrawableVectorField V;
    DrawableVectorField Vgrid;
    DrawableVectorField GT;

    std::vector<double> areas;

    GT=compute_ground_truth(grid,method,a,b,10,Eggs);
    F=get_scalar_field(m,a,b,Eggs,0,0);
    V=compute_field(m,F,method,2);


    bring_the_field_inside(m,grid,V,Vgrid,method);

    double err=estimate_MSE(GT,Vgrid,grid,method,without_boundary,mode);

    return err;
}

//
//===------------------------- EGGS FUNCTION ------------------------------------===//
//
//      f(x,y)=asin(bx)cos(by)
//
//===----------------------------------------------------------------------===//
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
int main(int argc, char *argv[])
{

    QApplication a(argc, argv);
    std::vector<double> Np(2);
    std::vector<double> anis(3);
    std::vector<double> Na(3);
    std::vector<double> Nb(3);

    //Build the grid
    DrawableTrimesh<> grid;
    make_grid(grid,500);
    //To choose the sweep
    int boundary;
    int mode;
    int method;
    int triangulation;
    int noise=0;
    int sweep;

    //Fixed values(when needed)
    double anisotropy=6;
    double A=0.1;
    double B=10;
    int N=25;
    double k=1;



    /*std::string s = "/Users/Dala/Documents/GitHub/gradient_fields/data/SFERA_"+nameNp0+".obj";*/
    DrawableTrimesh<> m;
    //make_triangulation(m,Delaunay,N,anisotropy);





    std::cout<<"Select the mode for this sweep:"<<std::endl;
    std::cout<<"0 range the number of samples"<<std::endl;
    std::cout<<"1 vary frequency"<<std::endl;
    std::cout<<"2 vary anisotropy"<<std::endl;
    std::cout<<"3 profiling"<<std::endl;
    std::cin>>mode;


    switch(mode)
    {
    case 0:
    {

        std::cout<<"Select the initial value of the frequency parameter:"<<std::endl;
        std::cin>>anis[0];
        std::cout<<"Select the final value of the frequency parameter:"<<std::endl;
        std::cin>>anis[1];
        std::cout<<"Select the number of samples for the frequency parameter:"<<std::endl;
        std::cin>>anis[2];

        assert(anis[1]>anis[0]);



        double delta_anis=(anis[1]-anis[0])/anis[2];
        std::ofstream outfile;

        std::ostringstream streamAnis0;
        streamAnis0 << std::fixed;
        streamAnis0 << std::setprecision(1);
        streamAnis0 << anis[0];
        std::string anis0=streamAnis0.str();

        std::ostringstream streamAnis1;
        streamAnis1 << std::fixed;
        streamAnis1 << std::setprecision(1);
        streamAnis1 << anis[1];
        std::string anis1=streamAnis1.str();

        outfile.open("/Users/Dala/Documents/MATLAB/Gradient_fields/TRI/PCEvsWorld.csv");
        for(int i=0;i<anis[2]+1;++i)
        {
            make_triangulation(m,0,anis[0]+delta_anis*i,anisotropy);



            for(int j=5;j<36;++j)
            {
                std::vector<double> values(4);
                for(int s=0;s<4;++s)
                {

                    ScalarField F=get_scalar_field(m,A,j,Eggs,0,0);
                    double err=compute(m,grid,s,A,j,3,N,F,anisotropy,1,1,0);
                    values[s]=err;
                    //std::vector<double> err=compute(m,grid,h,A,B,1,N,anis[0]+delta_anis*i,1,total,0,1);



                    //                    for (int j=0;j<err.size();++j)
                    //                    {
                    //                        if(err[j]==0)
                    //                        {continue;}
                    //                        else{outfile<<err[j]<<",";}
                    //                    }
                    //                     outfile<<'\n';


                }
                if(values[0]<values[1] && values[0]<values[2] && values[0]<values[3])
                {
                    outfile<<j<<","<<std::endl;
                    break;

                }else{
                    continue;
                }


            }
        }
        std::cout<<"FINISH"<<std::endl;
    }



        break;
    case 1:
    {


        std::cout<<"Select the initial value of the frequency parameter:"<<std::endl;
        std::cin>>anis[0];
        std::cout<<"Select the final value of the frequency parameter:"<<std::endl;
        std::cin>>anis[1];
        std::cout<<"Select the number of samples for the frequency parameter:"<<std::endl;
        std::cin>>anis[2];

        assert(anis[1]>anis[0]);



        std::cout<<"Select the sweep:"<<std::endl;

        std::cout<<"0 All the methods"<<std::endl;
        std::cout<<"1 Choose the method"<<std::endl;
        std::cin>>sweep;
        if(sweep)
        {
            std::cout<<"Select the method:"<<std::endl;
            std::cin>>method;

            //            std::cout<<"Which type of vertices do you want to consider?"<<std::endl;
            //            std::cout<<"0 All of them"<<std::endl;
            //            std::cout<<"1 Exclude boundary vertices "<<std::endl;
            //            std::cout<<"2 Consider ONLY boundary vertices"<<std::endl;

            //            std::cin>>boundary;


            double delta_anis=(anis[1]-anis[0])/anis[2];

            for(int v=0;v<2;++v)
            {
                N=(v==0)?25:50;
                make_triangulation(m,v,N,anisotropy);

                for (int h=0;h<3;++h)

                {
                    std::ofstream outfile;

                    std::ostringstream streamAnis0;
                    streamAnis0 << std::fixed;
                    streamAnis0 << std::setprecision(1);
                    streamAnis0 << anis[0];
                    std::string anis0=streamAnis0.str();

                    std::ostringstream streamAnis1;
                    streamAnis1 << std::fixed;
                    streamAnis1 << std::setprecision(1);
                    streamAnis1 << anis[1];
                    std::string anis1=streamAnis1.str();



                    outfile.open("/Users/Dala/Documents/MATLAB/Gradient_fields/TRI/VaryingFreq/Eggs_new/relative/"+ anis0 + "to"+anis1+" method"+method_names[method] +mode_names[h]+tri_names[v]+"noise.csv");

                    for(int i=0;i<anis[2]+1;++i)
                    {
                        //std::vector<double> err=compute(m,grid,method,A,B,3,N,anis[0]+delta_anis*i,1,total,0,1);
                        ScalarField F;//=get_scalar_field(m,A,B,Eggs,1,anis[0]+delta_anis*i);
                        double err=compute(m,grid,method,A,anis[0]+delta_anis*i,3,N,F,anisotropy,1,h,2);
                        outfile<<err<<","<<std::endl;
                    }

                }

            std::cout<<"FINISH"<<std::endl;
        }else
        {
            /*std::cout<<"Which type of vertices do you want to consider?"<<std::endl;
            std::cout<<"0 All of them"<<std::endl;
            std::cout<<"1 Exclude boundary vertices "<<std::endl;
            std::cout<<"2 Consider ONLY boundary vertices"<<std::endl;

            std::cin>>boundary;*/


            double delta_anis=(anis[1]-anis[0])/anis[2];
            for(int s=0;s<4;++s)
            {


                for(int h=0;h<1;++h)
                {


                    std::ofstream outfile;

                    std::ostringstream streamAnis0;
                    streamAnis0 << std::fixed;
                    streamAnis0 << std::setprecision(1);
                    streamAnis0 << anis[0];
                    std::string anis0=streamAnis0.str();



                    std::ostringstream streamAnis1;
                    streamAnis1 << std::fixed;
                    streamAnis1 << std::setprecision(1);
                    streamAnis1 << anis[1];
                    std::string anis1=streamAnis1.str();

                    outfile.open("/Users/Dala/Documents/MATLAB/Gradient_fields/TRI/VaryingFreq/relative/"+ anis0 + "to"+anis1+" method"+method_names[s]+mode_names[h]+"Boundary.csv");


                    for(int i=0;i<anis[2]+1;++i)
                    {


                        ScalarField F;//=get_scalar_field(m,A,B,Eggs,1,anis[0]+delta_anis*i);
                        double err=compute(m,grid,s,A,anis[0]+delta_anis*i,3,N,F,anisotropy,2,total,2);


                        //std::vector<double> err=compute(m,grid,h,A,B,1,N,anis[0]+delta_anis*i,1,total,0,1);



                        //                    for (int j=0;j<err.size();++j)
                        //                    {
                        //                        if(err[j]==0)
                        //                        {continue;}
                        //                        else{outfile<<err[j]<<",";}
                        //                    }
                        //                     outfile<<'\n';

                        outfile<<err<<","<<std::endl;
                    }

                }
                std::cout<<"FINISH with"+method_names[s]+"now"<<std::endl;
            }
            std::cout<<"FINISH"<<std::endl;

        }
    }
        break;
    case 2:
    {


        std::vector<double> frequencies={1,5,9};


        for(int k=0;k<3;++k)
        {
            std::ofstream outfile;
            std::ostringstream streamB;
            streamB << std::fixed;
            streamB << std::setprecision(1);
            streamB << frequencies[k];
            std::string B_s=streamB.str();

            outfile.open("/Users/Dala/Documents/MATLAB/Gradient_fields/TRI/box_plots/Distribution"+B_s+"anis.csv");
            make_triangulation(m,3,N,frequencies[k]);
            for(int s=0;s<4;++s)
            {
                ScalarField F;//=get_scalar_field(m,A,B,Eggs,1,frequencies[k]);

                std::vector<double> err=compute_for_hist(m,grid,s,A,B,1,N,F,anisotropy,1,total,2);



                for (int j=0;j<err.size();++j)
                {

                    outfile<<err[j]<<",";
                }
                outfile<<'\n';

            }
            std::cout<<"FINISH with B="+B_s+"now"<<std::endl;
        }
        std::cout<<"FINISH"<<std::endl;

    }
        break;
    case 3:
    {
        for(int method=0;method<2;++method)
        {
            std::ofstream outfile;
            outfile.open("/Users/Dala/Documents/MATLAB/Gradient_fields/ProfilingOnSphereof"+ method_names[method] +" .csv");

            for(int i=0;i<6;++i)
            {
                std::ostringstream streamNp0;
                streamNp0 << std::fixed;
                streamNp0 << std::setprecision(0);
                streamNp0 << i;
                std::string nameNp0=streamNp0.str();
                std::string s = "/Users/Dala/Documents/GitHub/gradient_fields/data/QUAD_"+nameNp0+".off";
                DrawableTrimesh<> m(s.c_str());

                std::cout<<m.num_verts()<<std::endl;
                ScalarField f=get_scalar_field(m,0.1,10,Eggs,0,1);
                //if (method!=0 && method !=1)
                //{
                std::vector<double> areas(m.num_polys());
                for(int k=0;k<m.num_polys();++k)
                {
                    areas[k]=m.poly_area(k);
                }

                std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

                DrawableVectorField V=compute_field(m,f,method);

                std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);

                outfile<<time_span.count()<<","<<std::endl;
                /*}else if(method==0)
                {
                    std::chrono::high_resolution_clock::time_point t0matrix = std::chrono::high_resolution_clock::now();
                    Eigen::SparseMatrix<double> G=gradient_matrix(m,1);
                    std::chrono::high_resolution_clock::time_point t1matrix = std::chrono::high_resolution_clock::now();
                    std::chrono::duration<double> time_spanmatrix = std::chrono::duration_cast<std::chrono::duration<double>>(t1matrix - t0matrix);

                    outfile<<time_spanmatrix.count()<<",";


                    /*DrawableVectorField V=DrawableVectorField(m,1);
                    std::chrono::high_resolution_clock::time_point t0mult = std::chrono::high_resolution_clock::now();
                    V=G*f;
                    std::chrono::high_resolution_clock::time_point t1mult = std::chrono::high_resolution_clock::now();
                    std::chrono::duration<double> time_spanmult = std::chrono::duration_cast<std::chrono::duration<double>>(t1mult - t0mult);

                    outfile<<time_spanmult.count()<<","<<std::endl;
                }else
                {
                    std::chrono::high_resolution_clock::time_point t0matrixags = std::chrono::high_resolution_clock::now();
                    Eigen::SparseMatrix<double> G=gradient_matrix(m,1);
                    std::chrono::high_resolution_clock::time_point t1matrixags = std::chrono::high_resolution_clock::now();
                    std::chrono::duration<double> time_spanmatrixags = std::chrono::duration_cast<std::chrono::duration<double>>(t1matrixags - t0matrixags);

                    outfile<<time_spanmatrixags.count()<<",";

                    DrawableVectorField W=DrawableVectorField(m,1);

                    DrawableVectorField Wv=DrawableVectorField(m,0);
                    vec3d provv=vec3d(0,0,0);
                    vec3d avg=vec3d(0,0,0);
                    float a=0;

                    std::vector<uint> nbr;

                   std::chrono::high_resolution_clock::time_point t0mult = std::chrono::high_resolution_clock::now();
                   W=G*f;
                    for (int k=0;k<m.num_verts();++k)
                    {
                        nbr=m.adj_v2p(k);
                        provv=vec3d(0,0,0);
                        avg=vec3d(0,0,0);
                        a=0;
                        for(int j=0;j<nbr.size();++j)
                        {
                            uint pid=nbr[j];
                            provv=W.vec_at(pid).operator *(m.poly_area(pid));
                            a+=m.poly_area(pid);
                            avg=avg.operator +(provv);

                        }
                        avg=avg.operator /(a);
                        Wv.set(k,avg);
                    }
                    std::chrono::high_resolution_clock::time_point t1mult = std::chrono::high_resolution_clock::now();
                    std::chrono::duration<double> time_spanmult = std::chrono::duration_cast<std::chrono::duration<double>>(t1mult - t0mult);

                    /*DrawableVectorField V=DrawableVectorField(m,0);
                    std::chrono::high_resolution_clock::time_point t0mult = std::chrono::high_resolution_clock::now();
                    V=G*f;
                    std::chrono::high_resolution_clock::time_point t1mult = std::chrono::high_resolution_clock::now();
                    std::chrono::duration<double> time_spanmult = std::chrono::duration_cast<std::chrono::duration<double>>(t1mult - t0mult);*/

                // outfile<<time_spanmult.count()<<","<<std::endl;
            }

        }

    }
        std::cout<<"FINISH"<<std::endl;







        break;





    }

    return a.exec();
}
