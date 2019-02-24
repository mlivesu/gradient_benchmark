#include <QApplication>
#include <cinolib/meshes/meshes.h>
#include <cinolib/gui/qt/qt_gui_tools.h>
#include "gui.h"
#include<cinolib/profiler.h>
using namespace cinolib;

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/*ScalarField compute( DrawableTetmesh<> &m,  DrawableTetmesh<> &grid,const int method,const double a, const double b,const int tet_type, const int N,const double anisotropy,int without_boundary,int mode,int noise,double k)
{


    make_tet_mesh(m,N,tet_type);
    DrawableVectorField V;
    DrawableVectorField Vgrid;
    DrawableVectorField GT;


    if(method!=0)
    {
        V=DrawableVectorField(m,false);
        GT=DrawableVectorField(m,false);

    }else
    {
        V=DrawableVectorField(m,true);
        GT=DrawableVectorField(m,true);
    }

    ScalarField F=get_scalar_field(m,a,b,Paraboloid,noise,k);
    GT=compute_ground_truth(m,method,a,b,10,Paraboloid);
    V=compute_field(m,F,method);
    //bring_the_field_inside(m,grid,V,Vgrid,method);
    ScalarField err=estimate_error(GT,V,m,method,without_boundary,mode);
    return err;
}*/
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
COMPUTATIONS computations_cubic;
int main(int argc, char **argv)
{
    QApplication a(argc, argv);


    GUI gui;
    init_gui(gui);
    init_events(gui);


    make_tet_mesh(computations_cubic.m,13,0);

//
    make_tet_mesh(computations_cubic.m_grid,50,0);
    gui.canvas.push_obj(&computations_cubic.m);
//    profile.pop();





       std::cout<<computations_cubic.m.num_verts()<<std::endl;
       std::cout<<computations_cubic.m.edge_avg_length()<<std::endl;

//       std::cout<<computations_cubic.m_grid.num_verts()<<std::endl;
//       std::cout<<computations_cubic.m_grid.edge_avg_length()<<std::endl;




//    // CMD+1 to show mesh controls.
//    VolumeMeshControlPanel<DrawableTetmesh<>> panel(&computations_cubic.m, &gui);
//    QApplication::connect(new QShortcut(QKeySequence(Qt::CTRL+Qt::Key_1), &gui), &QShortcut::activated, [&](){panel.show();});

    return a.exec();
}
