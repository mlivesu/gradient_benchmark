#include <QApplication>
#include<iostream>
#include<iomanip>
#include <cinolib/meshes/meshes.h>
#include<cinolib/profiler.h>
#include <cinolib/gui/qt/qt_gui_tools.h>
#include "meshes.h"
#include "computations_cubic.h"
#include "functions_cubic.h"

using namespace cinolib;
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double compute( DrawableTetmesh<> &m,  DrawableTetmesh<> &grid,const int method,const double a, const double b,const double c,const int tet_type, const int N,ScalarField &F,const double anisotropy,int without_boundary,int mode,int noise,int k)

{



    DrawableVectorField V;
    DrawableVectorField Vgrid=DrawableVectorField(grid,false);
    DrawableVectorField GT;

    GT=compute_ground_truth(grid,method,a,b,10,Eggs);
    F=get_scalar_field(m,a,b,10,Eggs,noise,k);
    V=compute_field(m,F,method,2);
    bring_the_field_inside(m,grid,V,Vgrid,method);
    double err=estimate_MSE(GT,Vgrid,grid,method,without_boundary,mode);

    return err;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
int main(int argc, char *argv[])
{
    QApplication a(argc, argv);

    //Parameter to vary

    std::vector<double> anis(3);


    //MAKE GRID
    DrawableTetmesh<> grid;
    make_tet_mesh(grid,100,0);



    //To choose the sweep
    int boundary;
    int mode;
    int method;
    int sweep;
    int mesh;
    int tetrahelization;

    //Fixed values(when needed)
    double anisotropy=6;
    double A=0.1;
    double B=10;
    double C=10;
    int N=50;
    double k=1;
    std::chrono::duration<double> time_precom;
    std::chrono::duration<double> time_estimation;


    DrawableTetmesh<> m;
    ScalarField F;
    Eigen::SparseMatrix<double> G;
    DrawableVectorField V;





    std::cout<<"Select the mode for this sweep:"<<std::endl;
    std::cout<<"0 varying frequency"<<std::endl;
    std::cout<<"1 varying anisotropy" <<std::endl;
    std::cout<<"2 varying noise"<<std::endl;
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

        std::cout<<"Select the mesh:"<<std::endl;

        std::cout<<"0 Regular"<<std::endl;
        std::cout<<"1 Delaunay"<<std::endl;

        std::cin>>mesh;

        std::cout<<"Which type of vertices do you want to consider?"<<std::endl;
        std::cout<<"0 All of them"<<std::endl;
        std::cout<<"1 Exclude boundary vertices "<<std::endl;
        std::cout<<"2 Consider ONLY boundary vertices"<<std::endl;

        std::cin>>boundary;

        std::cout<<"Select the sweep:"<<std::endl;

        std::cout<<"0 All the methods"<<std::endl;
        std::cout<<"1 Choose the method"<<std::endl;
        std::cin>>sweep;
        if(sweep)
        {
            std::cout<<"Select the method:"<<std::endl;
            std::cin>>method;

            N=(mesh==0)?25:50;

            make_tet_mesh(m,mesh,N,anisotropy);
            double delta_anis=(anis[1]-anis[0])/anis[2];

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



                outfile.open("Type the path of the folder in which you want to save the result of your experiment");

                for(int i=0;i<anis[2]+1;++i)
                {


                    double err=compute(m,grid,method,A,anis[0]+delta_anis*i,C,3,N,F,anisotropy,boundary,h,0,1);
                    outfile<<err<<","<<std::endl;
                }

            }

            std::cout<<"FINISH"<<std::endl;
        }else
        {


            double delta_anis=(anis[1]-anis[0])/anis[2];
            N=(mesh==0)?25:50;

            make_tet_mesh(m,mesh,N,anisotropy);

            for(int s=1;s<4;++s)
            {


                for(int h=0;h<3;++h)
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

                    outfile.open("Type the path of the folder in which you want to save the result of your experiment");


                    for(int i=0;i<anis[2]+1;++i)
                    {
                        double err=compute(m,grid,method,A,anis[0]+delta_anis*i,C,3,N,F,anisotropy,boundary,h,0,1);
                        outfile<<err<<","<<std::endl;
                    }

                }
                std::cout<<"FINISH with"+method_names[s]+"now"<<std::endl;
            }

            std::cout<<"FINISH"<<std::endl;

        }
    }



        break;
    case 1:
    {


        std::cout<<"Select the initial value of the anisotropy parameter:"<<std::endl;
        std::cin>>anis[0];
        std::cout<<"Select the final value of the anisotropy parameter:"<<std::endl;
        std::cin>>anis[1];
        std::cout<<"Select the number of samples for the anisotropy parameter:"<<std::endl;
        std::cin>>anis[2];

        assert(anis[1]>anis[0]);


        std::cout<<"Which type of vertices do you want to consider?"<<std::endl;
        std::cout<<"0 All of them"<<std::endl;
        std::cout<<"1 Exclude boundary vertices "<<std::endl;
        std::cout<<"2 Consider ONLY boundary vertices"<<std::endl;

        std::cin>>boundary;

        std::cout<<"Select the sweep:"<<std::endl;

        std::cout<<"0 All the methods"<<std::endl;
        std::cout<<"1 Choose the method"<<std::endl;
        std::cin>>sweep;
        if(sweep)
        {
            std::cout<<"Select the method:"<<std::endl;
            std::cin>>method;




            double delta_anis=(anis[1]-anis[0])/anis[2];

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



                outfile.open("Type the path of the folder in which you want to save the result of your experiment");

                for(int i=0;i<anis[2]+1;++i)
                {

                    make_tet_mesh(m,2,N,anis[0]+delta_anis*i);
                    double err=compute(m,grid,method,A,B,C,3,N,F,anisotropy,boundary,h,0,1);
                    outfile<<err<<","<<std::endl;
                }

            }
            std::cout<<"FINISH"<<std::endl;
        }
        else
        {


            double delta_anis=(anis[1]-anis[0])/anis[2];



            for(int s=1;s<4;++s)
            {


                for(int h=0;h<3;++h)
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

                    outfile.open("Type the path of the folder in which you want to save the result of your experiment");


                    for(int i=0;i<anis[2]+1;++i)
                    {
                        make_tet_mesh(m,2,N,anis[0]+delta_anis*i);
                        double err=compute(m,grid,method,A,B,C,3,N,F,anisotropy,boundary,h,0,1);
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

        std::cout<<"Select the initial value of the noise parameter:"<<std::endl;
        std::cin>>anis[0];
        std::cout<<"Select the final value of the noise parameter:"<<std::endl;
        std::cin>>anis[1];
        std::cout<<"Select the number of samples for the noise parameter:"<<std::endl;
        std::cin>>anis[2];

        assert(anis[1]>anis[0]);

        std::cout<<"Select the mesh:"<<std::endl;

        std::cout<<"0 Regular"<<std::endl;
        std::cout<<"1 Delaunay"<<std::endl;

        std::cin>>mesh;

        std::cout<<"Which type of vertices do you want to consider?"<<std::endl;
        std::cout<<"0 All of them"<<std::endl;
        std::cout<<"1 Exclude boundary vertices "<<std::endl;
        std::cout<<"2 Consider ONLY boundary vertices"<<std::endl;

        std::cin>>boundary;

        std::cout<<"Select the sweep:"<<std::endl;

        std::cout<<"0 All the methods"<<std::endl;
        std::cout<<"1 Choose the method"<<std::endl;
        std::cin>>sweep;
        if(sweep)
        {
            std::cout<<"Select the method:"<<std::endl;
            std::cin>>method;

            N=(mesh==0)?25:50;

            make_tet_mesh(m,mesh,N,anisotropy);

            double delta_anis=(anis[1]-anis[0])/anis[2];

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



                outfile.open("Type the path of the folder in which you want to save the result of your experiment");

                for(int i=0;i<anis[2]+1;++i)
                {


                    double err=compute(m,grid,method,A,B,C,3,N,F,anisotropy,boundary,h,1,anis[0]+delta_anis*i);
                    outfile<<err<<","<<std::endl;
                }

            }
            std::cout<<"FINISH"<<std::endl;
        }
        else
        {


            double delta_anis=(anis[1]-anis[0])/anis[2];
            make_triangulation(m,mesh,N,anisotropy);


            for(int s=0;s<4;++s)
            {


                for(int h=0;h<3;++h)
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

                    outfile.open("Type the path of the folder in which you want to save the result of your experiment");


                    for(int i=0;i<anis[2]+1;++i)
                    {

                        double err=compute(m,grid,method,A,B,C,3,N,F,anisotropy,boundary,h,1,anis[0]+delta_anis*i);
                        outfile<<err<<","<<std::endl;
                    }

                }
                std::cout<<"FINISH with"+method_names[s]+"now"<<std::endl;
            }

            std::cout<<"FINISH"<<std::endl;

        }


    }
        break;
    case 3:
    {

        for(int method=0;method<4;++method)
        {
            std::ofstream outfile;
            outfile.open("Type the path of the folder in which you want to save the result of your experiment");

            for(int i=0;i<6;++i)
            {
                std::ostringstream streamNp0;
                streamNp0 << std::fixed;
                streamNp0 << std::setprecision(0);
                streamNp0 << i;
                std::string nameNp0=streamNp0.str();
                std::string s; //=Type the path to the meshes on which you want to test the methods
                DrawableTrimesh<> m(s.c_str());

                switch(method)
                {
                case 0:
                {
                    V=DrawableVectorField(m,true);
                    ScalarField f=get_scalar_field(m,0.1,10,10,Eggs,0,1);
                    std::chrono::high_resolution_clock::time_point t0mult = std::chrono::high_resolution_clock::now();
                    G=gradient_matrix(m);
                    std::chrono::high_resolution_clock::time_point t1mult = std::chrono::high_resolution_clock::now();

                    time_precom = std::chrono::duration_cast<std::chrono::duration<double>>(t1mult - t0mult);


                    std::chrono::high_resolution_clock::time_point t0sol = std::chrono::high_resolution_clock::now();
                    V=G*f;
                    std::chrono::high_resolution_clock::time_point t1sol = std::chrono::high_resolution_clock::now();
                    time_estimation = std::chrono::duration_cast<std::chrono::duration<double>>(t1sol - t0sol);

                    outfile<<time_precom.count()<<",";

                    outfile<<time_estimation.count()<<","<<std::endl;
                }break;

                case 1:
                {
                    V=DrawableVectorField(m,false);
                    ScalarField f=get_scalar_field(m,0.1,10,10,Eggs,0,1);
                    std::chrono::high_resolution_clock::time_point t0mult = std::chrono::high_resolution_clock::now();
                    G=build_matrix_for_AGS(m);
                    std::chrono::high_resolution_clock::time_point t1mult = std::chrono::high_resolution_clock::now();

                    time_precom = std::chrono::duration_cast<std::chrono::duration<double>>(t1mult - t0mult);


                    std::chrono::high_resolution_clock::time_point t0sol = std::chrono::high_resolution_clock::now();
                    V=G*f;
                    std::chrono::high_resolution_clock::time_point t1sol = std::chrono::high_resolution_clock::now();
                    time_estimation = std::chrono::duration_cast<std::chrono::duration<double>>(t1sol - t0sol);

                    outfile<<time_precom.count()<<",";

                    outfile<<time_estimation.count()<<","<<std::endl;
                }
                    break;
                case 2:
                {

                    std::vector<Eigen::ColPivHouseholderQR<Eigen::MatrixXd> > M;
                    std::vector<Eigen::MatrixXd> RHS;
                    std::vector<std::vector<uint>> nbrs;

                    ScalarField f=get_scalar_field(m,0.1,10,10,Eggs,0,1);

                    std::chrono::high_resolution_clock::time_point t0mult = std::chrono::high_resolution_clock::now();
                    build_matrix_for_LSDD(m,M,RHS,nbrs);
                    std::chrono::high_resolution_clock::time_point t1mult = std::chrono::high_resolution_clock::now();

                    time_precom = std::chrono::duration_cast<std::chrono::duration<double>>(t1mult - t0mult);


                    std::chrono::high_resolution_clock::time_point t0sol = std::chrono::high_resolution_clock::now();
                    solve_for_LSDD(m,V,M,RHS,f,nbrs);
                    std::chrono::high_resolution_clock::time_point t1sol = std::chrono::high_resolution_clock::now();
                    time_estimation = std::chrono::duration_cast<std::chrono::duration<double>>(t1sol - t0sol);

                    outfile<<time_precom.count()<<",";

                    outfile<<time_estimation.count()<<","<<std::endl;
                }
                    break;
                case 3:
                {
                    DrawableVectorField V;
                    std::vector<Eigen::ColPivHouseholderQR<Eigen::MatrixXd> > M;
                    std::vector<Eigen::MatrixXd> RHS;
                    std::vector<std::vector<uint>> nbrs;
                    ScalarField f=get_scalar_field(m,0.1,10,10,Eggs,0,1);

                    std::chrono::high_resolution_clock::time_point t0mult = std::chrono::high_resolution_clock::now();

                    build_matrix_for_LR(m,M,RHS,nbrs);

                    std::chrono::high_resolution_clock::time_point t1mult = std::chrono::high_resolution_clock::now();
                    time_precom = std::chrono::duration_cast<std::chrono::duration<double>>(t1mult - t0mult);

                    std::chrono::high_resolution_clock::time_point t0sol = std::chrono::high_resolution_clock::now();
                    solve_for_LR(m,V,M,RHS,f,nbrs);

                    std::chrono::high_resolution_clock::time_point t1sol = std::chrono::high_resolution_clock::now();
                    time_estimation = std::chrono::duration_cast<std::chrono::duration<double>>(t1sol - t0sol);

                    outfile<<time_precom.count()<<",";

                    outfile<<time_estimation.count()<<","<<std::endl;
                }

                break;

            }
        }
        std::cout<<"FINISH"<<std::endl;


    }
    }
        break;
    }

    return a.exec();
}
