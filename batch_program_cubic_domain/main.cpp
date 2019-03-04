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
double compute( DrawableTetmesh<> &m,  DrawableTetmesh<> &grid,const int method,const double a, const double b,const int tet_type, const int N,ScalarField &F,const double anisotropy,int without_boundary,int mode,int weight)
//std::vector<double> compute( DrawableTetmesh<> &m,  DrawableTetmesh<> &grid,const int method,const double a, const double b,const int tet_type, const int N,const double anisotropy,int without_boundary,int mode,int noise,double k)
{


    //make_tet_mesh(m,N,3,anisotropy);
    DrawableVectorField V;
    DrawableVectorField Vgrid=DrawableVectorField(grid,false);
    DrawableVectorField GT;

    GT=compute_ground_truth(grid,method,a,b,10,Eggs);
    //F=get_scalar_field(m,a,b,10,Eggs,0,0);
    V=compute_field(m,F,method,2);
    bring_the_field_inside(m,grid,V,Vgrid,method);
    double err=estimate_MSE(GT,Vgrid,grid,method,without_boundary,mode);
    //std::vector<double> err=error_for_hist(GT,Vgrid,grid,method,without_boundary,mode);
    return err;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
int main(int argc, char *argv[])
{
    QApplication a(argc, argv);

    //RANGES FOR THE PARAMETERS
    std::vector<double> Np(2);
    std::vector<double> anis(3);
    std::vector<double> Na(3);
    std::vector<double> Nb(3);


    //MAKE GRID
    DrawableTetmesh<> grid;
    //make_tet_mesh(grid,100,0);



    //To choose the sweep
    int boundary;
    int mode;
    int method;
    int sweep;
    int tetrahelization;

    //Fixed values(when needed)
    double anisotropy=6;
    double A=0.1;
    double B=10;
    int N=50;
    double k=1;
    std::chrono::duration<double> time_precom;
    std::chrono::duration<double> time_estimation;


//    //MAKE MESH(WHEN FIXED)
//    DrawableTetmesh<> m;

//    //make_tet_mesh(m,N,1);



//    std::cout<<"Select the mode for this sweep:"<<std::endl;
//    std::cout<<"0 range the number of samples"<<std::endl;
//    std::cout<<"1 vary frequency "<<std::endl;
//    std::cout<<"2 vary anisotropy"<<std::endl;
//    std::cout<<"3 profiling"<<std::endl;
//    std::cin>>mode;

//    switch(mode)
//    {
//    case 0:
//    {
//        std::cout<<"Select the initial number of samples you want to start with:"<<std::endl;
//        std::cin>>Np[0];
//        assert(!std::cin.fail());
//        std::cout<<"Select the final number of samples you want to end with:"<<std::endl;
//        std::cin>>Np[1];
//        assert(!std::cin.fail());
//        assert(Np[1]>Np[0]);

//        std::cout<<"Select the method:"<<std::endl;
//        std::cin>>method;
//        std::ofstream outfile;

//        std::ostringstream streamNp0;
//        streamNp0 << std::fixed;
//        streamNp0 << std::setprecision(1);
//        streamNp0 << Np[0];
//        std::string nameNp0=streamNp0.str();

//        std::ostringstream streamNp1;
//        streamNp1 << std::fixed;
//        streamNp1 << std::setprecision(1);
//        streamNp1 << Np[1];
//        std::string nameNp1=streamNp1.str();
//        outfile.open("/Users/Dala/Documents/MATLAB/Gradient_fields/Nfrom"+ nameNp0 + "to" + nameNp1 +" .csv");
//        /*  make_tet_mesh(m,N,1);
//        std::vector<uint>nbr;
//        double tmp;
//        double count=0;

//        for(int i=0;i<m.num_verts();++i)
//        {
//            if(m.vert_is_on_srf(i))
//            {
//                ++count;
//                continue;
//            }
//            else{
//                 tmp+=m.vert_valence(i);
//            }


//        }

//        std::cout<<tmp/(m.num_verts()-count)<<std::endl;*/

//        for(int i=0;i<Np[1]-Np[0];++i)
//        {
//            //double err=compute(m,grid,method,A,B,1,Np[0]+i,anisotropy,1,total,0,1);
//            //            std::vector<double> err=compute(m,grid,method,A,B,1,Np[0]+i,anisotropy,1,total,0,1);
//            //            for (int j=0;j<err.size();++j)
//            //            {
//            //                outfile<<err[j]<<",";
//            //            }
//            //             outfile<<std::endl;
//            //            //outfile<<err<<","<<std::endl;
//        }
//        std::cout<<"FINISH"<<std::endl;
//    }
//        break;
//    case 1:
//    {
//        std::cout<<"Select the initial value for the frequency:"<<std::endl;
//        std::cin>>anis[0];
//        std::cout<<"Select the final value for the frequency:"<<std::endl;
//        std::cin>>anis[1];
//        std::cout<<"Select the number of samples for the frequency:"<<std::endl;
//        std::cin>>anis[2];

//        assert(anis[1]>anis[0]);

//        std::cout<<"Select the sweep:"<<std::endl;

//        std::cout<<"0 All the methods"<<std::endl;
//        std::cout<<"1 Choose the method"<<std::endl;
//        std::cin>>sweep;
//        if(sweep)
//        {
//            std::cout<<"Select the method:"<<std::endl;
//            std::cin>>method;
//            double delta_anis=(anis[1]-anis[0])/anis[2];
////            for(int v=0;v<2;++v)

////            {
////                N=(v==0)?25:50;
////                make_tet_mesh(m,N,v);
//                for (int h=1;h<3;++h)

//                {
//                    std::ofstream outfile;

//                    std::ostringstream streamAnis0;
//                    streamAnis0 << std::fixed;
//                    streamAnis0 << std::setprecision(1);
//                    streamAnis0 << anis[0];
//                    std::string anis0=streamAnis0.str();

//                    std::ostringstream streamAnis1;
//                    streamAnis1 << std::fixed;
//                    streamAnis1 << std::setprecision(1);
//                    streamAnis1 << anis[1];
//                    std::string anis1=streamAnis1.str();



//                    outfile.open("/Users/Dala/Documents/MATLAB/Gradient_fields/TET/IncrementingNoise/relative/"+ anis0 + "to"+anis1+" method"+method_names[method] +mode_names[h]+".csv");

//                    for(int i=0;i<anis[2]+1;++i)
//                    {
//                        //std::vector<double> err=compute(m,grid,method,A,B,3,N,anis[0]+delta_anis*i,1,total,0,1);
//                        ScalarField F=get_scalar_field(m,A,B,10,Eggs,1,anis[0]+delta_anis*i);
//                        double err=compute(m,grid,method,A,B,3,N,F,anisotropy,1,h,2);
//                        outfile<<err<<","<<std::endl;
//                    }
//                    std::cout<<"FINISH"<<std::endl;
//                }
//            //}
//        }else
//        {
//            double delta_anis=(anis[1]-anis[0])/anis[2];

////            for(int v=1;v<2;++v)
////            {
////                N=(v==0)?25:50;
////                make_tet_mesh(m,N,v);

//            for(int s=2;s<3;++s)
//            {


//                for(int h=1;h<3;++h)
//                {

//                    std::ofstream outfile;
//                    std::ofstream outfile_hist;

//                    std::ostringstream streamAnis0;
//                    streamAnis0 << std::fixed;
//                    streamAnis0 << std::setprecision(1);
//                    streamAnis0 << anis[0];
//                    std::string anis0=streamAnis0.str();

//                    std::ostringstream streamB;
//                    streamB << std::fixed;
//                    streamB << std::setprecision(1);
//                    streamB << B;
//                    std::string B_s=streamB.str();

//                    std::ostringstream streamAnis1;
//                    streamAnis1 << std::fixed;
//                    streamAnis1 << std::setprecision(1);
//                    streamAnis1 << anis[1];
//                    std::string anis1=streamAnis1.str();

//                    outfile.open("/Users/Dala/Documents/MATLAB/Gradient_fields/TET/IncrementingNoise/relative/"+ anis0 + "to"+anis1+" method"+method_names[s]+mode_names[h]+".csv");



//                    for(int i=0;i<anis[2]+1;++i)
//                    {
//                        ScalarField F=get_scalar_field(m,A,B,10,Eggs,1,anis[0]+delta_anis*i);
//                        double err=compute(m,grid,s,A,B,3,N,F,anisotropy,1,h,2);
//                        //std::vector<double> err=compute(m,grid,h,A,B,1,N,anis[0]+delta_anis*i,1,total,0,1);



//                        //                    for (int j=0;j<err.size();++j)
//                        //                    {
//                        //                        if(err[j]==0)
//                        //                        {continue;}
//                        //                        else{outfile<<err[j]<<",";}
//                        //                    }
//                        //                     outfile<<'\n';

//                        outfile<<err<<","<<std::endl;
//                    }
//                    std::cout<<"FINISH with"+method_names[s]+mode_names[h]+"now"<<std::endl;
//                }
//            }
//            //}
//            std::cout<<"FINISH"<<std::endl;
//        }





//    }break;
//    case 2:
//    {

//        std::vector<double> frequencies={1,5,9};
//        double scale=1.87;
//        //      for(int v=0;v<2;++v)

//        //      {
//        //          N=(v==0)?13:25;
//        //make_tet_mesh(m,25,1);

//        for(int k=0;k<3;++k)
//        {
//            std::ofstream outfile;
//            std::ostringstream streamB;
//            streamB << std::fixed;
//            streamB << std::setprecision(1);
//            streamB << frequencies[k];
//            std::string B_s=streamB.str();

//            outfile.open("/Users/Dala/Documents/MATLAB/Gradient_fields/TET/box_plots/DistributionAnis"+B_s+".csv");
//            make_tet_mesh(m,13,3,frequencies[k]);
//            for(int s=0;s<4;++s)
//            {
//                ScalarField F;//=get_scalar_field(m,A*scale,10/scale,10,Eggs,1,frequencies[k]);

//                std::vector<double> err=compute_for_hist(m,grid,s,A*scale,10/scale,1,N,F,frequencies[k],1,total,2);



//                for (int j=0;j<err.size();++j)
//                {

//                    outfile<<err[j]<<",";
//                }
//                outfile<<'\n';

//            }
//            std::cout<<"FINISH with B="+ B_s +" now"<<std::endl;
//        }
//        //}
//        std::cout<<"FINISH"<<std::endl;


//    }break;
//    case 3:
//    {
//        Profiler profiler;

        for(int method=3;method<4;++method)
        {
            std::ofstream outfile;
            outfile.open("/Users/Dala/Documents/MATLAB/Gradient_fields/Profiling/3D/Onlyone"+ method_names[method] +" .csv");

            for(int i=5;i<6;++i)
            {

                std::ostringstream streamNp0;
                streamNp0 << std::fixed;
                streamNp0 << std::setprecision(0);
                streamNp0 << i;
                std::string nameNp0=streamNp0.str();
                std::string s = "/Users/Dala/Documents/GitHub/gradient_fields/data/CUBO_"+nameNp0+".MESH";
                DrawableTetmesh<> m(s.c_str());


                if(method==2)
                {
                    DrawableVectorField V;
                    std::vector<Eigen::ColPivHouseholderQR<Eigen::MatrixXd> > M;
                    std::vector<Eigen::MatrixXd> RHS;
                    std::vector<std::vector<uint>> nbrs;

                    ScalarField f=get_scalar_field(m,A,B,10,Eggs,1,0);

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
                }else
                {
                    DrawableVectorField V;
                    std::vector<Eigen::ColPivHouseholderQR<Eigen::MatrixXd> > M;
                    std::vector<Eigen::MatrixXd> RHS;
                    std::vector<std::vector<uint>> nbrs;
                    ScalarField f=get_scalar_field(m,A,B,10,Eggs,1,0);

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



        }
        }
        std::cout<<"FINISH"<<std::endl;

//    }

//    }
//         break;
//    }

    return a.exec();
}
