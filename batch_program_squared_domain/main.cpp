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


     //make_triangulation(m,tri_type,N,anisotropy);

    DrawableVectorField V;
    DrawableVectorField Vgrid;
    DrawableVectorField GT;

    std::vector<double> areas;

    GT=compute_ground_truth(grid,method,a,b,10,Eggs);
    //F=get_scalar_field(m,a,b,Eggs,0,0);
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
    //make_grid(grid,500);
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
    std::chrono::duration<double> time_precom;
    std::chrono::duration<double> time_estimation;



    DrawableTrimesh<> m;
    //make_triangulation(m,Delaunay,N,anisotropy);



    std::cout<<"Select the mode for this sweep:"<<std::endl;
    std::cout<<"0 valence experiment"<<std::endl;
    std::cout<<"1 vary frequency"<<std::endl;
    std::cout<<"2 vary anisotropy"<<std::endl;
    std::cout<<"3 profiling"<<std::endl;
    std::cin>>mode;


    switch(mode)
    {
    case 0:
    {


        DrawableVectorField V;
        DrawableVectorField Vgrid;
        DrawableVectorField GT;
        for (int s=1;s<4;++s)
        {
            std::ofstream outfile;

            std::vector<std::vector<double>> matrix(31, std::vector<double>(10));




            for(int b=5;b<36;++b)
            {


                GT=compute_ground_truth(m,s,A,b,10,Eggs);
                ScalarField F=get_scalar_field(m,A,b,Eggs,0,0);
                V=compute_field(m,F,s,2);
                std::vector<double> errors(10);
                std::vector<int> count(10);
                for(int i=0;i<m.num_verts();++i)
                {
                    if(!m.vert_is_boundary(i))
                    {
                        int valence=m.vert_valence(i);
                        double err=relative_error(GT.vec_at(i),V.vec_at(i),total);

                        errors[valence]=errors[valence]+err;
                        count[valence]=count[valence]+1;
                    }


                }
                for(uint k=0;k<errors.size();++k)
                {

                    matrix[b-5][k]=errors[k]/count[k];
                    std::cout<<k<<std::endl;
                    std::cout<<count[k]<<std::endl;

                }
            }

            outfile.open("/Users/Dala/Documents/MATLAB/Gradient_fields/TRI/valence/"+ method_names[s]+".csv");

            for(uint v=0;v<31;++v)
            {
                for(uint w=0;w<10;++w)
                {
                    outfile<<matrix[v][w]<<",";
                }
                outfile<<'\n';
            }
        }

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



                    outfile.open("/Users/Dala/Documents/MATLAB/Gradient_fields/TRI/IncrementingNoise/relative/"+ anis0 + "to"+anis1+" method"+method_names[method] +mode_names[h]+tri_names[v]+"noise.csv");

                    for(int i=0;i<anis[2]+1;++i)
                    {
                        //std::vector<double> err=compute(m,grid,method,A,B,3,N,anis[0]+delta_anis*i,1,total,0,1);
                        ScalarField F;//=get_scalar_field(m,A,B,Eggs,1,anis[0]+delta_anis*i);
                        double err=compute(m,grid,method,A,anis[0]+delta_anis*i,3,N,F,anisotropy,1,h,2);
                        outfile<<err<<","<<std::endl;
                    }

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

//            for(int v=0;v<2;++v)
//            {

//                make_triangulation(m,v,25,anisotropy);
                for(int s=2;s<4;++s)
                {


                    for(int h=1;h<3;++h)
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

                        outfile.open("/Users/Dala/Documents/MATLAB/Gradient_fields/TRI/IncrementingNoise/relative/"+ anis0 + "to"+anis1+" method"+method_names[s]+mode_names[h]+".csv");


                        for(int i=0;i<anis[2]+1;++i)
                        {


                            ScalarField F=get_scalar_field(m,A,B,Eggs,1,anis[0]+delta_anis*i);
                            double err=compute(m,grid,s,A,B,3,N,F,anisotropy,1,h,2);


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
                //std::cout<<"FINISH with"+tri_names[v]+"now"<<std::endl;
            //}
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

        for(int method=2;method<3;++method)
        {
            std::ofstream outfile;
            outfile.open("/Users/Dala/Documents/MATLAB/Gradient_fields/Profiling/2D/ProfilingOnSphereof"+ method_names[method] +" .csv");

            for(int i=0;i<6;++i)
            {
                std::ostringstream streamNp0;
                streamNp0 << std::fixed;
                streamNp0 << std::setprecision(0);
                streamNp0 << i;
                std::string nameNp0=streamNp0.str();
                std::string s = "/Users/Dala/Documents/GitHub/gradient_fields/data/QUAD_"+nameNp0+".off";
                DrawableTrimesh<> m(s.c_str());


                if(method==2)
                {
                    DrawableVectorField V;
                    std::vector<Eigen::ColPivHouseholderQR<Eigen::Matrix2d> > M;
                    std::vector<Eigen::MatrixXd> RHS;
                    std::vector<std::vector<uint>> nbrs;

                    ScalarField f=get_scalar_field(m,0.1,10,Eggs,0,1);

                    std::chrono::high_resolution_clock::time_point t0mult = std::chrono::high_resolution_clock::now();
                    build_matrix_for_LSDD(m,M,RHS,nbrs);
                    std::chrono::high_resolution_clock::time_point t1mult = std::chrono::high_resolution_clock::now();

                    time_precom = std::chrono::duration_cast<std::chrono::duration<double>>(t1mult - t0mult);


                    std::chrono::high_resolution_clock::time_point t0sol = std::chrono::high_resolution_clock::now();
                    solve_for_LSDD(m,V,M,RHS,f,nbrs,time_precom,time_estimation);
                    std::chrono::high_resolution_clock::time_point t1sol = std::chrono::high_resolution_clock::now();
                    time_estimation = std::chrono::duration_cast<std::chrono::duration<double>>(t1sol - t0sol);

                    outfile<<time_precom.count()<<",";

                    outfile<<time_estimation.count()<<","<<std::endl;
                }else
                {
                    DrawableVectorField V;
                    std::vector<Eigen::ColPivHouseholderQR<Eigen::MatrixXd> > M;
                    std::vector<Eigen::MatrixXd> RHS;
                    std::vector<std::vector<uint>> nbrs;
                    ScalarField f=get_scalar_field(m,0.1,10,Eggs,0,1);

                    std::chrono::high_resolution_clock::time_point t0mult = std::chrono::high_resolution_clock::now();

                    build_matrix_for_LR(m,M,RHS,nbrs);

                    std::chrono::high_resolution_clock::time_point t1mult = std::chrono::high_resolution_clock::now();
                    time_precom = std::chrono::duration_cast<std::chrono::duration<double>>(t1mult - t0mult);

                    std::chrono::high_resolution_clock::time_point t0sol = std::chrono::high_resolution_clock::now();
                    solve_for_LR(m,V,M,RHS,f,nbrs,time_precom,time_estimation);

                    std::chrono::high_resolution_clock::time_point t1sol = std::chrono::high_resolution_clock::now();
                    time_estimation = std::chrono::duration_cast<std::chrono::duration<double>>(t1sol - t0sol);

                    outfile<<time_precom.count()<<",";

                    outfile<<time_estimation.count()<<","<<std::endl;
                }

            }

        }
        std::cout<<"FINISH"<<std::endl;


    }




        break;





    }

    return a.exec();
}
