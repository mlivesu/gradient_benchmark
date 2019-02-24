#include <QGridLayout>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QLabel>
#include "gui.h"
#include "functions.h"
#include <QDir>
#include<iostream>
#include<fstream>
#include "picking.h"


extern COMPUTATIONS computations;

using namespace cinolib;

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void init_gui(GUI & gui)
{
    QGridLayout *global_layout  = new QGridLayout();
    QGridLayout *choose_heat_map_layout=new QGridLayout();
    QGroupBox *show= new QGroupBox;
    QGridLayout *error_layout=new QGridLayout();
    QGridLayout *show_layout=new QGridLayout();

    //CANVAS
    gui.canvas.set_2d_mode(true);

    //PANEL
    QGridLayout *button_layout= new QGridLayout();
    QVBoxLayout *toolbar_layout = new QVBoxLayout();


    //PANEL CANVAS I
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    button_layout->addWidget(new QLabel("Analytic Function:"),0,0);
    button_layout->addWidget(&gui.cb_function,0,1);
    gui.cb_function.setMaximumSize(200,50);
    button_layout->addWidget(new QLabel("10*a"),1,0);
    button_layout->addWidget(&gui.a,1,1);



    gui.a.setMaximumSize(50,50);

    gui.a.setRange(1,50);
    gui.a.setValue(1);

    button_layout->addWidget(new QLabel("b"),1,2);
    button_layout->addWidget(&gui.b,1,3);

    gui.b.setMaximumSize(100,50);
    gui.b.setRange(-50,50);
    gui.b.setValue(10);

    button_layout->addWidget(new QLabel("c"),1,4);
    button_layout->addWidget(&gui.c,1,5);
    gui.c.setMaximumSize(100,50);
    gui.c.setRange(-50,50);
    gui.c.setValue(10);

    button_layout->addWidget(new QLabel("Methods:"),2,0);
    button_layout->addWidget(&gui.method,2,1);
    gui.method.setMaximumSize(200,50);
    button_layout->addWidget(new QLabel("Triangulation:"),3,0);
    button_layout->addWidget(&gui.choose_tri,3,1);
    gui.choose_tri.setMaximumSize(150,50);
    button_layout->addWidget(new QLabel("Number of samples:"),4,0);
    button_layout->addWidget(&gui.N,4,1);
    gui.N.setMaximumSize(100,50);
    gui.N.setMinimum(1);
    gui.N.setMaximum(500);
    gui.N.setValue(20);

    button_layout->addWidget(new QLabel("Anisotropy:"),5,0);
    button_layout->addWidget(&gui.anis,5,1);
    gui.anis.setMaximumSize(100,50);
    gui.anis.setMinimum(1);
    gui.anis.setMaximum(9);
    gui.anis.setValue(4);


    show_layout->addWidget(&gui.but,0,0);
    show_layout->addWidget(&gui.ground_truth,0,1);
    show_layout->addWidget(&gui.show_heat_map,1,0);
    show_layout->addWidget(&gui.wireframe,1,1);
    gui.wireframe.setChecked(true);
    show->setTitle("Show");
    show->setLayout(show_layout);

    button_layout->addWidget(show,6,1);

    choose_heat_map_layout->addWidget(&gui.scalar_field,0,0);
    choose_heat_map_layout->addWidget(&gui.Err,0,1);

    gui.choose_heat_map.setTitle("Choose Heat Map");
    gui.choose_heat_map.setLayout(choose_heat_map_layout);
    gui.choose_heat_map.setDisabled(true);

    button_layout->addWidget(&gui.choose_heat_map,9,1);

    error_layout->addWidget(new QLabel("Mode:"),0,0);
    error_layout->addWidget(&gui.mode,0,1);
    error_layout->addWidget(new QLabel("Error type:"),1,0);
    error_layout->addWidget(&gui.ErrType,1,1);
    error_layout->addWidget(new QLabel("Interpolated Error:"),2,0);
    error_layout->addWidget(&gui.InsideError,2,1);
    error_layout->addWidget(new QLabel("Type of Vertices:"),3,0);
    error_layout->addWidget(&gui.type_of_vertices,3,1);
    error_layout->addWidget(new QLabel("Handle Boundaries:"),4,0);
    error_layout->addWidget(&gui.boundaries,4,1);
    gui.error.setTitle("Error");
    gui.error.setLayout(error_layout);
    gui.error.setDisabled(true);

    button_layout->addWidget(&gui.error,13,1);






    button_layout->addWidget(&gui.save,14,0);
    button_layout->addWidget(&gui.reset,14,1);
//    button_layout->addWidget(new QLabel("Rescale Scalar Field"),10,0);
//    gui.sl_scalar_field.setOrientation(Qt::Horizontal);

//    button_layout->addWidget(&gui.sl_scalar_field,10,1);
//    gui.sl_scalar_field.setRange(1,50);
//    gui.sl_scalar_field.setValue(25);


    button_layout->addWidget(new QLabel("Negative Sautration"),10,0);
    button_layout->addWidget(new QLabel("Positive Sautration"),11,0);
    gui.sl_error_neg.setOrientation(Qt::Horizontal);
    gui.sl_error_pos.setOrientation(Qt::Horizontal);

    button_layout->addWidget(&gui.sl_error_neg,10,1);
    button_layout->addWidget(&gui.sl_error_pos,11,1);
    gui.sl_error_neg.setRange(0,500);
    gui.sl_error_neg.setValue(0);

    gui.sl_error_pos.setRange(250,1000);
    gui.sl_error_pos.setValue(1000);

    /*button_layout->addWidget(new QLabel("Rescale Vector Field"),9,2);
    button_layout->addWidget(&gui.sl_vector_fields,9,3);

    gui.sl_vector_fields.setRange(1,10);
    gui.sl_vector_fields.setValue(1);
    gui.sl_vector_fields.setDisabled(true);*/

    button_layout->setSpacing(0);
    button_layout->setHorizontalSpacing(0);
    button_layout->setContentsMargins(0,0,0,0);
    toolbar_layout->addLayout(button_layout);
    toolbar_layout->setSpacing(0);

    toolbar_layout->setContentsMargins(0,0,0,0);


    //

    for(int i=0; i<N_FUNCTIONS; ++i) gui.cb_function.addItem(f_names[i].c_str());
    for(int i=0; i<=LR_CENTROIDS; ++i) gui.method.addItem(method_names[i].c_str());
    for(int i=0; i<=grid; ++i) gui.choose_tri.addItem(tri_names[i].c_str());
    for(int i=0; i<=relative; ++i) gui.ErrType.addItem(err_names[i].c_str());
    for(int i=0; i<=angle; ++i) gui.mode.addItem(mode_names[i].c_str());
    for(int i=0; i<=only_boundary; ++i) gui.type_of_vertices.addItem(vertices_names[i].c_str());
    gui.but.setText("Gradient");
    gui.ground_truth.setText("Ground Truth");
    gui.wireframe.setText("Wireframe");
    gui.scalar_field.setText("Scalar Field");
    gui.Err.setText("Error");
    gui.save.setText("Save");
    gui.reset.setText("Reset");
    gui.show_heat_map.setText("Show Heat Map");

    toolbar_layout->addStretch();


    global_layout->addLayout(toolbar_layout,0,0);


    gui.widget.setLayout(global_layout);
    gui.widget.show();
    gui.canvas.show();
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void init_events(GUI & gui)
{
    gui.canvas.setMinimumSize(500,500);
    // bind mouse click handler
    gui.canvas.callback_mouse_press = [&](GLcanvas *c, QMouseEvent *e)
    {
        if(e->modifiers() == Qt::ControlModifier)
        {

            vec3d p;
            vec2i click(e->x(), e->y());
            if (c->unproject(click, p)) // transform click in a 3d point
            {
                uint vid=closest_vertex(p,&computations.m);
                Eigen::VectorXd X=fit_with_quadrics(computations.m,vid,computations.f);
                for(int i=0;i<6;++i)
                {
                    std::cout<<X(i)<<std::endl;
                }
//                if(gui.InsideError.isChecked())
//                {
//                    vid=closest_vertex(p,&computations.m_grid);
//                    std::cout << "vid in the grid:" << std::endl;
//                    std::cout << vid << std::endl;
//                }else if (gui.Err.isChecked())
//                {
//                    if(gui.method.currentIndex()!=0)
//                    {
//                        vid= closest_vertex(p,&computations.m);
//                    }else
//                    {
//                        vid = closest_vertex(p,&computations.dual_m);
//                    }
//                    std::cout << mode_names[gui.mode.currentIndex()].c_str() <<" error calculated:"<<std::endl;
//                    std::cout<<computations.err[vid] << std::endl;
//                }else
//                {
//                    vid= closest_vertex(p,&computations.m);
//                    std::cout << computations.f[vid] << std::endl;
//                    //std::cout << vid << std::endl;
//                }
            }
        }
        else if(e->modifiers() == Qt::ShiftModifier)
        {

            vec3d p;
            vec2i click(e->x(), e->y());
            if (c->unproject(click, p)) // transform click in a 3d point
            {

                uint vid = closest_vertex(p,&computations.m);
                std::cout<<"Function value="<<computations.f[vid]<<std::endl;
                std::cout<<"Vert position="<<std::endl;
                std::cout<<computations.m.vert(vid)[0]<<std::endl;
                std::cout<<computations.m.vert(vid)[1]<<std::endl;
                std::cout<<computations.m.vert(vid)[2]<<std::endl;
            }
            c->updateGL();


        }
    };


    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    QCheckBox::connect(&gui.show_heat_map, &QCheckBox::stateChanged, [&]()
    {
        if(gui.show_heat_map.isChecked())
        {

            gui.choose_heat_map.setDisabled(false);
            gui.sl_scalar_field.setDisabled(false);
            gui.sl_error_neg.setDisabled(false);
            gui.sl_error_pos.setDisabled(false);

            if (gui.scalar_field.isChecked())
            {
                gui.canvas.pop(&computations.dual_m);
                gui.canvas.pop(&computations.m_grid);

                computations.f=get_scalar_field(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());
                find_max_min_values(computations.f,computations.max,computations.min);
                computations.f_norm=heat_map_normalization(computations.f,computations.min,computations.max);
                computations.f_norm.copy_to_mesh(computations.m);
                computations.m.show_texture1D(TEXTURE_1D_HSV);


            }

            if (gui.Err.isChecked())
            {
                gui.error.setDisabled(false);

                gui.canvas.pop(&computations.dual_m);
                gui.canvas.pop(&computations.m_grid);

                if(gui.InsideError.isChecked())
                {
                    computations.f=get_scalar_field(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());
                    computations.V=compute_field(computations.m,computations.f,gui.method.currentIndex());
                    computations.GT_grid=compute_values_on_grid(computations.m_grid, gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());

                    bring_the_field_inside(computations.m,computations.m_grid,computations.V,computations.V_grid,gui.method.currentIndex());
                    computations.GT=compute_ground_truth(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());

                    computations.err=estimate_error(computations.GT_grid,computations.V_grid,computations.m_grid,gui.mode.currentIndex(),LSDD,gui.ErrType.currentIndex(),gui.type_of_vertices.currentIndex());
                    computations.max=find_max_norm(computations.GT);
                    computations.err_norm=heat_map_normalization(computations.err,0,computations.max,gui.sl_error_neg.value(),gui.sl_error_pos.value());
                    computations.err_norm.copy_to_mesh(computations.m_grid);

                    gui.canvas.pop(&computations.m);
                    gui.canvas.push_obj(&computations.m_grid,false);
                    computations.m_grid.show_wireframe(false);
                    computations.m_grid.show_texture1D(TEXTURE_1D_HSV);
                    gui.canvas.push_obj(&computations.m,false);

                }else
                {
                    computations.f=get_scalar_field(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());
                    computations.GT=compute_ground_truth(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
                    computations.V=compute_field(computations.m,computations.f,gui.method.currentIndex());
                    computations.err=estimate_error(computations.GT,computations.V,computations.m,gui.mode.currentIndex(),gui.method.currentIndex(),gui.ErrType.currentIndex(),gui.type_of_vertices.currentIndex());
                    computations.max=find_max_norm(computations.GT);
                    computations.err_norm=heat_map_normalization(computations.err,0,computations.max,gui.sl_error_neg.value(),gui.sl_error_pos.value());

                    if(gui.method.currentIndex()==0)
                    {
                        computations.err_norm.copy_to_mesh(computations.dual_m);
                        gui.canvas.pop(&computations.m);
                        gui.canvas.push_obj(&computations.dual_m,false);
                        computations.dual_m.show_wireframe(false);
                        computations.dual_m.show_texture1D(TEXTURE_1D_HSV);
                        gui.canvas.push_obj(&computations.m,false);

                    }else
                    {
                        computations.err_norm.copy_to_mesh(computations.m);
                        computations.m.show_texture1D(TEXTURE_1D_HSV);


                    }
                }



                if(gui.but.isChecked())
                {
                    computations.f=get_scalar_field(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());
                    computations.V=compute_field(computations.m,computations.f,gui.method.currentIndex());
                    computations.V_norm=computations.V;
                    computations.V_norm.normalize();
                    gui.canvas.push_obj(&computations.V_norm,false);


                }
                if(gui.ground_truth.isChecked())
                {
                    computations.GT=compute_ground_truth(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
                    computations.GT_norm=computations.GT;
                    computations.GT_norm.normalize();
                    gui.canvas.push_obj(&computations.GT_norm,false);

                }

            }

        }else{

            gui.choose_heat_map.setDisabled(true);
            gui.sl_scalar_field.setDisabled(true);
            gui.sl_error_neg.setDisabled(true);
            gui.sl_error_pos.setDisabled(true);
            gui.error.setDisabled(true);
            computations.dual_m.show_vert_color();
            computations.m_grid.show_vert_color();
            computations.m.show_vert_color();
        }


        gui.canvas.updateGL();
    });
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    QCheckBox::connect(&gui.but, &QCheckBox::stateChanged, [&]()
    {
        if(gui.but.isChecked())
        {
            computations.f=get_scalar_field(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());
            computations.V=compute_field(computations.m,computations.f,gui.method.currentIndex());
            computations.V_norm=computations.V;
            computations.V_norm.normalize();
            gui.canvas.push_obj(&computations.V_norm,false);


        }
        else
        {
            gui.canvas.pop_all_occurrences_of(DRAWABLE_VECTOR_FIELD);

            if(gui.ground_truth.isChecked())
            {
                computations.GT=compute_ground_truth(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
                computations.GT_norm=computations.GT;
                computations.GT_norm.normalize();
                gui.canvas.push_obj(&computations.GT_norm,false);

            }
        }
        gui.canvas.updateGL();



    });
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    QSlider::connect(&gui.sl_scalar_field, &QSlider::valueChanged, [&]()
    {
        if (gui.scalar_field.isChecked())
        {
                        gui.canvas.pop(&computations.dual_m);
                        gui.canvas.pop(&computations.m_grid);
            //            //computations.dual_m.clear();


                       computations.f=get_scalar_field(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());
                     find_max_min_values(computations.f,computations.max,computations.min);
            computations.f_norm=heat_map_normalization(computations.f,computations.min,computations.max);
            computations.f_norm.copy_to_mesh(computations.m);
            computations.m.show_texture1D(TEXTURE_1D_HSV);


        }

        if (gui.Err.isChecked())
        {
            /*gui.canvas.pop(&computations.dual_m);
            gui.canvas.pop(&computations.m_grid);*/
            if(gui.InsideError.isChecked())
            {
                computations.f=get_scalar_field(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());
                computations.GT=compute_ground_truth(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
                computations.V=compute_field(computations.m,computations.f,gui.method.currentIndex());
                computations.err=estimate_error(computations.GT,computations.V,computations.m,gui.mode.currentIndex(),gui.method.currentIndex(),gui.ErrType.currentIndex(),gui.type_of_vertices.currentIndex());
                computations.max=find_max_norm(computations.GT);
                computations.err_norm=heat_map_normalization(computations.err,0,computations.max,gui.sl_error_neg.value(),gui.sl_error_pos.value());
                computations.err_norm.copy_to_mesh(computations.m_grid);

                gui.canvas.pop(&computations.m);
                gui.canvas.push_obj(&computations.m_grid,false);
                computations.m_grid.show_wireframe(false);

                computations.m_grid.show_texture1D(TEXTURE_1D_HSV);
                computations.m.show_vert_color();

                gui.canvas.push_obj(&computations.m,false);

            }else{

                computations.f=get_scalar_field(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());
                computations.GT=compute_ground_truth(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
                computations.V=compute_field(computations.m,computations.f,gui.method.currentIndex());
                computations.err=estimate_error(computations.GT,computations.V,computations.m,gui.mode.currentIndex(),gui.method.currentIndex(),gui.ErrType.currentIndex(),gui.type_of_vertices.currentIndex());
                computations.max=find_max_norm(computations.GT);
                computations.err_norm=heat_map_normalization(computations.err,0,computations.max,gui.sl_error_neg.value(),gui.sl_error_pos.value());

                if(gui.method.currentIndex()==0)
                {


                    computations.err_norm.copy_to_mesh(computations.dual_m);

                    gui.canvas.pop(&computations.m);
                    gui.canvas.push_obj(&computations.dual_m,false);
                    computations.dual_m.show_wireframe(false);

                    computations.dual_m.show_texture1D(TEXTURE_1D_HSV);
                    gui.canvas.push_obj(&computations.m,false);
                }else
                {
                    computations.err_norm.copy_to_mesh(computations.m);
                    computations.m.show_texture1D(TEXTURE_1D_HSV);


                }
            }
        }


        gui.canvas.updateGL();

    });
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    QCheckBox::connect(&gui.boundaries, &QCheckBox::stateChanged, [&]()
    {
        if(gui.boundaries.isChecked())
        {
            if(gui.method.currentIndex()!=1)
            {
                computations.f=get_scalar_field(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());
                computations.GT=compute_ground_truth(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
                computations.V=compute_field(computations.m,computations.f,gui.method.currentIndex(),true);
            }

        else
        {
            std::map<uint, int> boundaries;
            computations.f=scalar_field_with_boundaries(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),boundaries);
            computations.GT=compute_ground_truth(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
            computations.V=compute_field(computations.m,computations.f,gui.method.currentIndex(),true);
        }
        }else
        {
            computations.f=get_scalar_field(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());
            computations.GT=compute_ground_truth(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
            computations.V=compute_field(computations.m,computations.f,gui.method.currentIndex(),false);
        }

        if(gui.but.isChecked())
        {
            computations.V_norm=computations.V;
            computations.V_norm.normalize();
            gui.canvas.push_obj(&computations.V_norm,false);
        }

        if(gui.ground_truth.isChecked())
        {
            computations.GT_norm=computations.GT;
            computations.GT_norm.normalize();
            gui.canvas.push_obj(&computations.GT_norm,false);

        }

        if (gui.show_heat_map.isChecked())
        {
            if (gui.scalar_field.isChecked())
            {
                gui.canvas.pop(&computations.dual_m);
                gui.canvas.pop(&computations.m_grid);
                find_max_min_values(computations.f,computations.max,computations.min);
                computations.f_norm=heat_map_normalization(computations.f,computations.min,computations.max);
                computations.f_norm.copy_to_mesh(computations.m);
                computations.m.show_texture1D(TEXTURE_1D_HSV);


            }

            if (gui.Err.isChecked())
            {

                gui.canvas.pop(&computations.dual_m);
                gui.canvas.pop(&computations.m_grid);
                computations.dual_m.clear();

                if(gui.InsideError.isChecked())
                {

                    computations.GT_grid=compute_values_on_grid(computations.m_grid, gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());

                    bring_the_field_inside(computations.m,computations.m_grid,computations.V,computations.V_grid,gui.method.currentIndex());
                    computations.GT=compute_ground_truth(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
                    computations.err=estimate_error(computations.GT_grid,computations.V_grid,computations.m_grid,gui.mode.currentIndex(),LSDD,gui.ErrType.currentIndex(),gui.type_of_vertices.currentIndex());
                    computations.max=find_max_norm(computations.GT);
                    computations.err_norm=heat_map_normalization(computations.err,0,computations.max,gui.sl_error_neg.value(),gui.sl_error_pos.value());
                    computations.err_norm.copy_to_mesh(computations.m_grid);

                    gui.canvas.pop(&computations.m);
                    gui.canvas.push_obj(&computations.m_grid,false);
                    computations.m_grid.show_wireframe(false);
                    computations.m_grid.show_texture1D(TEXTURE_1D_HSV);
                    gui.canvas.push_obj(&computations.m,false);


                }else
                {

                    computations.err=estimate_error(computations.GT,computations.V,computations.m,gui.mode.currentIndex(),gui.method.currentIndex(),gui.ErrType.currentIndex(),gui.type_of_vertices.currentIndex());
                    computations.max=find_max_norm(computations.GT);
                    computations.err_norm=heat_map_normalization(computations.err,0,computations.max,gui.sl_error_neg.value(),gui.sl_error_pos.value());

                    if(gui.method.currentIndex()==0)
                    {
                        computations.err_norm.copy_to_mesh(computations.dual_m);
                        gui.canvas.pop(&computations.m);
                        gui.canvas.push_obj(&computations.dual_m,false);
                        computations.dual_m.show_wireframe(false);
                        computations.dual_m.show_texture1D(TEXTURE_1D_HSV);
                        gui.canvas.push_obj(&computations.m,false);


                    }else
                    {
                        computations.err_norm.copy_to_mesh(computations.m);
                        computations.m.show_texture1D(TEXTURE_1D_HSV);

                    }

                }



            }
        }

    });
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    QSlider::connect(&gui.sl_error_neg, &QSlider::valueChanged, [&]()
    {


        if (gui.Err.isChecked())
        {
            gui.canvas.pop(&computations.dual_m);
            gui.canvas.pop(&computations.m_grid);
            if(gui.InsideError.isChecked())
            {computations.f=get_scalar_field(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());
                computations.V=compute_field(computations.m,computations.f,gui.method.currentIndex());
                computations.GT_grid=compute_values_on_grid(computations.m_grid, gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());

                bring_the_field_inside(computations.m,computations.m_grid,computations.V,computations.V_grid,gui.method.currentIndex());
                computations.GT=compute_ground_truth(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());

                computations.err=estimate_error(computations.GT_grid,computations.V_grid,computations.m_grid,gui.mode.currentIndex(),LSDD,gui.ErrType.currentIndex(),gui.type_of_vertices.currentIndex());

                computations.max=find_max_norm(computations.GT);
                computations.err_norm=heat_map_normalization(computations.err,0,computations.max,gui.sl_error_neg.value(),gui.sl_error_pos.value());
                computations.err_norm.copy_to_mesh(computations.m_grid);

                gui.canvas.pop(&computations.m);
                gui.canvas.push_obj(&computations.m_grid,false);
                computations.m_grid.show_wireframe(false);
                computations.m_grid.show_texture1D(TEXTURE_1D_HSV);
                gui.canvas.push_obj(&computations.m,false);

            }else{

                computations.f=get_scalar_field(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());
                computations.GT=compute_ground_truth(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
                computations.V=compute_field(computations.m,computations.f,gui.method.currentIndex());
                computations.err=estimate_error(computations.GT,computations.V,computations.m,gui.mode.currentIndex(),gui.method.currentIndex(),gui.ErrType.currentIndex(),gui.type_of_vertices.currentIndex());
                computations.max=find_max_norm(computations.GT);
                computations.err_norm=heat_map_normalization(computations.err,0,computations.max,gui.sl_error_neg.value(),gui.sl_error_pos.value());

                if(gui.method.currentIndex()==0)
                {


                    computations.err_norm.copy_to_mesh(computations.dual_m);

                    gui.canvas.pop(&computations.m);
                    gui.canvas.push_obj(&computations.dual_m,false);
                    computations.dual_m.show_wireframe(false);


                    computations.dual_m.show_texture1D(TEXTURE_1D_HSV);
                    gui.canvas.push_obj(&computations.m,false);
                }else
                {
                    computations.err_norm.copy_to_mesh(computations.m);
                    computations.m.show_texture1D(TEXTURE_1D_HSV);


                }
            }
        }


        gui.canvas.updateGL();

    });
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    QSlider::connect(&gui.sl_error_pos, &QSlider::valueChanged, [&]()
    {


        if (gui.Err.isChecked())
        {
            gui.canvas.pop(&computations.dual_m);
            gui.canvas.pop(&computations.m_grid);
            if(gui.InsideError.isChecked())
            {
                computations.f=get_scalar_field(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());
                computations.V=compute_field(computations.m,computations.f,gui.method.currentIndex());
                computations.GT_grid=compute_values_on_grid(computations.m_grid, gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());

                bring_the_field_inside(computations.m,computations.m_grid,computations.V,computations.V_grid,gui.method.currentIndex());
                computations.GT=compute_ground_truth(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());

                computations.err=estimate_error(computations.GT_grid,computations.V_grid,computations.m_grid,gui.mode.currentIndex(),LSDD,gui.ErrType.currentIndex(),gui.type_of_vertices.currentIndex());
                computations.max=find_max_norm(computations.GT);
                computations.err_norm=heat_map_normalization(computations.err,0,computations.max,gui.sl_error_neg.value(),gui.sl_error_pos.value());
                computations.err_norm.copy_to_mesh(computations.m_grid);

                gui.canvas.pop(&computations.m);
                gui.canvas.push_obj(&computations.m_grid,false);
                computations.m_grid.show_wireframe(false);
                computations.m_grid.show_texture1D(TEXTURE_1D_HSV);
                gui.canvas.push_obj(&computations.m,false);

            }else{

                computations.f=get_scalar_field(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());
                computations.GT=compute_ground_truth(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
                computations.V=compute_field(computations.m,computations.f,gui.method.currentIndex());
                computations.err=estimate_error(computations.GT,computations.V,computations.m,gui.mode.currentIndex(),gui.method.currentIndex(),gui.ErrType.currentIndex(),gui.type_of_vertices.currentIndex());
                computations.max=find_max_norm(computations.GT);
                computations.err_norm=heat_map_normalization(computations.err,0,computations.max,gui.sl_error_neg.value(),gui.sl_error_pos.value());

                if(gui.method.currentIndex()==0)
                {


                    computations.err_norm.copy_to_mesh(computations.dual_m);

                    gui.canvas.pop(&computations.m);
                    gui.canvas.push_obj(&computations.dual_m,false);
                    computations.dual_m.show_wireframe(false);


                    computations.dual_m.show_texture1D(TEXTURE_1D_HSV);
                    gui.canvas.push_obj(&computations.m,false);
                }else
                {
                    computations.err_norm.copy_to_mesh(computations.m);
                    computations.m.show_texture1D(TEXTURE_1D_HSV);


                }
            }
        }


        gui.canvas.updateGL();

    });

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    QCheckBox::connect(&gui.ground_truth, &QCheckBox::stateChanged, [&]()
    {
        if(gui.ground_truth.isChecked())
        {
            computations.GT=compute_ground_truth(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
            computations.GT_norm=computations.GT;
            computations.GT_norm.normalize();
            gui.canvas.push_obj(&computations.GT_norm,false);

        }else
        {
            gui.canvas.pop_all_occurrences_of(DRAWABLE_VECTOR_FIELD);

            if(gui.but.isChecked())
            {
                computations.f=get_scalar_field(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());
                computations.V=compute_field(computations.m,computations.f,gui.method.currentIndex());
                computations.V_norm=computations.V;
                computations.V_norm.normalize();
                gui.canvas.push_obj(&computations.V_norm,false);


            }
        }


        gui.canvas.updateGL();
    });
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    QCheckBox::connect(&gui.wireframe, &QCheckBox::stateChanged, [&]()
    {
        if(gui.wireframe.isChecked())
        {
            computations.m.show_wireframe(true);
        }else{
            computations.m.show_wireframe(false);
        }
        gui.canvas.updateGL();
    });
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    QComboBox::connect(&gui.cb_function, static_cast<void(QComboBox::*)(int)>(&QComboBox::activated), [&]()
    {


        if(gui.but.isChecked())
        {
            computations.f=get_scalar_field(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());
            computations.V=compute_field(computations.m,computations.f,gui.method.currentIndex());
            computations.V_norm=computations.V;
            computations.V_norm.normalize();
            gui.canvas.push_obj(&computations.V_norm,false);


        }
        if(gui.ground_truth.isChecked())
        {
            computations.GT=compute_ground_truth(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
            computations.GT_norm=computations.GT;
            computations.GT_norm.normalize();
            gui.canvas.push_obj(&computations.GT_norm,false);

        }

        if (gui.show_heat_map.isChecked())
        {
            if (gui.scalar_field.isChecked())
            {
                gui.canvas.pop(&computations.dual_m);
                gui.canvas.pop(&computations.m_grid);



                computations.f=get_scalar_field(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());
                find_max_min_values(computations.f,computations.max,computations.min);
                computations.f_norm=heat_map_normalization(computations.f,computations.min,computations.max);
                computations.f_norm.copy_to_mesh(computations.m);
                computations.m.show_texture1D(TEXTURE_1D_HSV);


            }

            if (gui.Err.isChecked())
            {

                gui.canvas.pop(&computations.dual_m);
                gui.canvas.pop(&computations.m_grid);
                computations.dual_m.clear();

                if(gui.InsideError.isChecked())
                {

                    computations.f=get_scalar_field(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());
                    computations.V=compute_field(computations.m,computations.f,gui.method.currentIndex());
                    computations.GT_grid=compute_values_on_grid(computations.m_grid, gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());

                    bring_the_field_inside(computations.m,computations.m_grid,computations.V,computations.V_grid,gui.method.currentIndex());
                    computations.GT=compute_ground_truth(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
                    computations.err=estimate_error(computations.GT_grid,computations.V_grid,computations.m_grid,gui.mode.currentIndex(),LSDD,gui.ErrType.currentIndex(),gui.type_of_vertices.currentIndex());
                    computations.max=find_max_norm(computations.GT);
                    computations.err_norm=heat_map_normalization(computations.err,0,computations.max,gui.sl_error_neg.value(),gui.sl_error_pos.value());
                    computations.err_norm.copy_to_mesh(computations.m_grid);

                    gui.canvas.pop(&computations.m);
                    gui.canvas.push_obj(&computations.m_grid,false);
                    computations.m_grid.show_wireframe(false);
                    computations.m_grid.show_texture1D(TEXTURE_1D_HSV);
                    gui.canvas.push_obj(&computations.m,false);


                }else
                {
                    computations.f=get_scalar_field(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());
                    computations.GT=compute_ground_truth(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
                    computations.V=compute_field(computations.m,computations.f,gui.method.currentIndex());
                    computations.err=estimate_error(computations.GT,computations.V,computations.m,gui.mode.currentIndex(),gui.method.currentIndex(),gui.ErrType.currentIndex(),gui.type_of_vertices.currentIndex());
                    computations.max=find_max_norm(computations.GT);
                    computations.err_norm=heat_map_normalization(computations.err,0,computations.max,gui.sl_error_neg.value(),gui.sl_error_pos.value());

                    if(gui.method.currentIndex()==0)
                    {
                        computations.err_norm.copy_to_mesh(computations.dual_m);
                        gui.canvas.pop(&computations.m);
                        gui.canvas.push_obj(&computations.dual_m,false);
                        computations.dual_m.show_wireframe(false);
                        computations.dual_m.show_texture1D(TEXTURE_1D_HSV);
                        gui.canvas.push_obj(&computations.m,false);


                    }else
                    {
                        computations.err_norm.copy_to_mesh(computations.m);
                        computations.m.show_texture1D(TEXTURE_1D_HSV);

                    }

                }



            }
        }
        gui.canvas.updateGL();

    });
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    // defines reaction to a combo box item selection
    QComboBox::connect(&gui.method, static_cast<void(QComboBox::*)(int)>(&QComboBox::activated), [&]()
    {
        if(gui.but.isChecked())
        {
            computations.f=get_scalar_field(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());
            computations.V=compute_field(computations.m,computations.f,gui.method.currentIndex());
            computations.V_norm=computations.V;
            computations.V_norm.normalize();
            gui.canvas.push_obj(&computations.V_norm,false);


        }

        if(gui.ground_truth.isChecked())
        {
            computations.GT=compute_ground_truth(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
            computations.GT_norm=computations.GT;
            computations.GT_norm.normalize();
            gui.canvas.push_obj(&computations.GT_norm,false);

        }

        if (gui.show_heat_map.isChecked())
        {

            if (gui.scalar_field.isChecked())
            {
                gui.canvas.pop(&computations.dual_m);
                gui.canvas.pop(&computations.m_grid);



                computations.f=get_scalar_field(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());
                find_max_min_values(computations.f,computations.max,computations.min);
                computations.f_norm=heat_map_normalization(computations.f,computations.min,computations.max);
                computations.f_norm.copy_to_mesh(computations.m);
                computations.m.show_texture1D(TEXTURE_1D_HSV);


            }
            if (gui.Err.isChecked())
            {

                gui.canvas.pop(&computations.dual_m);
                gui.canvas.pop(&computations.m_grid);
                computations.dual_m.clear();

                if(gui.InsideError.isChecked())
                {

                    computations.f=get_scalar_field(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());
                    computations.V=compute_field(computations.m,computations.f,gui.method.currentIndex());
                    computations.GT_grid=compute_values_on_grid(computations.m_grid, gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());

                    bring_the_field_inside(computations.m,computations.m_grid,computations.V,computations.V_grid,gui.method.currentIndex());
                    computations.GT=compute_ground_truth(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
                    computations.err=estimate_error(computations.GT_grid,computations.V_grid,computations.m_grid,gui.mode.currentIndex(),LSDD,gui.ErrType.currentIndex(),gui.type_of_vertices.currentIndex());
                    computations.max=find_max_norm(computations.GT);
                    computations.err_norm=heat_map_normalization(computations.err,0,computations.max,gui.sl_error_neg.value(),gui.sl_error_pos.value());
                    computations.err_norm.copy_to_mesh(computations.m_grid);

                    gui.canvas.pop(&computations.m);
                    gui.canvas.push_obj(&computations.m_grid,false);
                    computations.m_grid.show_wireframe(false);
                    computations.m_grid.show_texture1D(TEXTURE_1D_HSV);
                    gui.canvas.push_obj(&computations.m,false);


                }else{

                    computations.f=get_scalar_field(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());
                    computations.GT=compute_ground_truth(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
                    computations.V=compute_field(computations.m,computations.f,gui.method.currentIndex());
                    computations.err=estimate_error(computations.GT,computations.V,computations.m,gui.mode.currentIndex(),gui.method.currentIndex(),gui.ErrType.currentIndex(),gui.type_of_vertices.currentIndex());
                    computations.max=find_max_norm(computations.GT);
                    computations.err_norm=heat_map_normalization(computations.err,0,computations.max,gui.sl_error_neg.value(),gui.sl_error_pos.value());

                    if(gui.method.currentIndex()==0)
                    {
                        computations.err_norm.copy_to_mesh(computations.dual_m);
                        gui.canvas.pop(&computations.m);
                        gui.canvas.push_obj(&computations.dual_m,false);
                        computations.dual_m.show_wireframe(false);
                        computations.dual_m.show_texture1D(TEXTURE_1D_HSV);
                        gui.canvas.push_obj(&computations.m,false);


                    }else
                    {
                        computations.err_norm.copy_to_mesh(computations.m);
                        computations.m.show_texture1D(TEXTURE_1D_HSV);

                    }
                }




            }
        }


        gui.canvas.updateGL();

    });
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    // defines reaction to a combo box item selection
    QComboBox::connect(&gui.choose_tri, static_cast<void(QComboBox::*)(int)>(&QComboBox::activated), [&]()
    {
        make_triangulation(computations.m,gui.choose_tri.currentIndex(),gui.N.value());


        if(gui.but.isChecked())
        {
            computations.f=get_scalar_field(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());
            computations.V=compute_field(computations.m,computations.f,gui.method.currentIndex());
            computations.V_norm=computations.V;
            computations.V_norm.normalize();
            gui.canvas.push_obj(&computations.V_norm,false);


        }
        if(gui.ground_truth.isChecked())
        {
            computations.GT=compute_ground_truth(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
            computations.GT_norm=computations.GT;
            computations.GT_norm.normalize();
            gui.canvas.push_obj(&computations.GT_norm,false);

        }

        if (gui.show_heat_map.isChecked())
        {

            if (gui.scalar_field.isChecked())
            {
                gui.canvas.pop(&computations.dual_m);
                gui.canvas.pop(&computations.m_grid);
                computations.dual_m.clear();


                computations.f=get_scalar_field(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());
                find_max_min_values(computations.f,computations.max,computations.min);
                computations.f_norm=heat_map_normalization(computations.f,computations.min,computations.max);
                computations.f_norm.copy_to_mesh(computations.m);
                computations.m.show_texture1D(TEXTURE_1D_HSV);


            }
            if (gui.Err.isChecked())
            {
                gui.canvas.pop(&computations.dual_m);
                gui.canvas.pop(&computations.m_grid);
                computations.dual_m.clear();
                if(gui.InsideError.isChecked())

               {
                    computations.f=get_scalar_field(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());
                    computations.V=compute_field(computations.m,computations.f,gui.method.currentIndex());
                    computations.GT_grid=compute_values_on_grid(computations.m_grid, gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());

                    bring_the_field_inside(computations.m,computations.m_grid,computations.V,computations.V_grid,gui.method.currentIndex());
                    computations.GT=compute_ground_truth(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
                    computations.err=estimate_error(computations.GT_grid,computations.V_grid,computations.m_grid,gui.mode.currentIndex(),LSDD,gui.ErrType.currentIndex(),gui.type_of_vertices.currentIndex());
                    computations.max=find_max_norm(computations.GT);
                    computations.err_norm=heat_map_normalization(computations.err,0,computations.max,gui.sl_error_neg.value(),gui.sl_error_pos.value());
                    computations.err_norm.copy_to_mesh(computations.m_grid);

                    gui.canvas.pop(&computations.m);
                    gui.canvas.push_obj(&computations.m_grid,false);
                    computations.m_grid.show_wireframe(false);
                    computations.m_grid.show_texture1D(TEXTURE_1D_HSV);
                    gui.canvas.push_obj(&computations.m,false);
                }
                else{




                    computations.f=get_scalar_field(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());
                    computations.GT=compute_ground_truth(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
                    computations.V=compute_field(computations.m,computations.f,gui.method.currentIndex());
                    computations.err=estimate_error(computations.GT,computations.V,computations.m,gui.mode.currentIndex(),gui.method.currentIndex(),gui.ErrType.currentIndex(),gui.type_of_vertices.currentIndex());
                    computations.max=find_max_norm(computations.GT);
                    computations.err_norm=heat_map_normalization(computations.err,0,computations.max,gui.sl_error_neg.value(),gui.sl_error_pos.value());

                    if(gui.method.currentIndex()==0)
                    {
                        computations.err_norm.copy_to_mesh(computations.dual_m);
                        gui.canvas.pop(&computations.m);
                        gui.canvas.push_obj(&computations.dual_m,false);
                        computations.dual_m.show_wireframe(false);
                        computations.dual_m.show_texture1D(TEXTURE_1D_HSV);
                        gui.canvas.push_obj(&computations.m,false);


                    }else
                    {
                        computations.err_norm.copy_to_mesh(computations.m);
                        computations.m.show_texture1D(TEXTURE_1D_HSV);

                    }
                }
            }
        }




        gui.canvas.updateGL();

    });
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    QSpinBox::connect(&gui.N,static_cast<void(QSpinBox::*)(int)>(&QSpinBox::valueChanged), [&]()
    {


        make_triangulation(computations.m,gui.choose_tri.currentIndex(),gui.N.value(),gui.anis.value());




        if(gui.but.isChecked())
        {
            computations.f=get_scalar_field(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());
            computations.V=compute_field(computations.m,computations.f,gui.method.currentIndex());
            computations.V_norm=computations.V;
            computations.V_norm.normalize();
            gui.canvas.push_obj(&computations.V_norm,false);


        }

        if(gui.ground_truth.isChecked())
        {
            computations.GT=compute_ground_truth(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
            computations.GT_norm=computations.GT;
            computations.GT_norm.normalize();
            gui.canvas.push_obj(&computations.GT_norm,false);

        }
        if (gui.show_heat_map.isChecked())
        {

            if (gui.scalar_field.isChecked())
            {


                computations.f=get_scalar_field(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());
                find_max_min_values(computations.f,computations.max,computations.min);
                computations.f_norm=heat_map_normalization(computations.f,computations.min,computations.max);
                computations.f_norm.copy_to_mesh(computations.m);
                computations.m.show_texture1D(TEXTURE_1D_HSV);


            }
            if (gui.Err.isChecked())
            {

                gui.canvas.pop(&computations.dual_m);
                gui.canvas.pop(&computations.m_grid);
                computations.dual_m.clear();

                if(gui.InsideError.isChecked())
                {

                    computations.f=get_scalar_field(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());
                    computations.V=compute_field(computations.m,computations.f,gui.method.currentIndex());
                    computations.GT_grid=compute_values_on_grid(computations.m_grid, gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());

                    bring_the_field_inside(computations.m,computations.m_grid,computations.V,computations.V_grid,gui.method.currentIndex());
                    computations.GT=compute_ground_truth(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
                    computations.err=estimate_error(computations.GT_grid,computations.V_grid,computations.m_grid,gui.mode.currentIndex(),LSDD,gui.ErrType.currentIndex(),gui.type_of_vertices.currentIndex());
                    computations.max=find_max_norm(computations.GT);
                    computations.err_norm=heat_map_normalization(computations.err,0,computations.max,gui.sl_error_neg.value(),gui.sl_error_pos.value());
                    computations.err_norm.copy_to_mesh(computations.m_grid);

                    gui.canvas.pop(&computations.m);
                    gui.canvas.push_obj(&computations.m_grid,false);
                    computations.m_grid.show_wireframe(false);
                    computations.m_grid.show_texture1D(TEXTURE_1D_HSV);
                    gui.canvas.push_obj(&computations.m,false);


                }else{


                    computations.f=get_scalar_field(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());
                    computations.GT=compute_ground_truth(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
                    computations.V=compute_field(computations.m,computations.f,gui.method.currentIndex());
                    computations.err=estimate_error(computations.GT,computations.V,computations.m,gui.mode.currentIndex(),gui.method.currentIndex(),gui.ErrType.currentIndex(),gui.type_of_vertices.currentIndex());
                    computations.max=find_max_norm(computations.GT);
                    computations.err_norm=heat_map_normalization(computations.err,0,computations.max,gui.sl_error_neg.value(),gui.sl_error_pos.value());

                    if(gui.method.currentIndex()==0)
                    {
                        computations.err_norm.copy_to_mesh(computations.dual_m);
                        gui.canvas.pop(&computations.m);
                        gui.canvas.push_obj(&computations.dual_m,false);
                        computations.dual_m.show_wireframe(false);
                        computations.dual_m.show_texture1D(TEXTURE_1D_HSV);
                        gui.canvas.push_obj(&computations.m,false);


                    }else
                    {
                        computations.err_norm.copy_to_mesh(computations.m);
                        computations.m.show_texture1D(TEXTURE_1D_HSV);

                    }
                }
            }
        }
        gui.canvas.updateGL();

    });
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    QSpinBox::connect(&gui.anis,static_cast<void(QSpinBox::*)(int)>(&QSpinBox::valueChanged), [&]()
    {
        make_triangulation(computations.m,gui.choose_tri.currentIndex(),gui.N.value(),gui.anis.value());

        if(gui.but.isChecked())
        {
            computations.f=get_scalar_field(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());
            computations.V=compute_field(computations.m,computations.f,gui.method.currentIndex());
            computations.V_norm=computations.V;
            computations.V_norm.normalize();
            gui.canvas.push_obj(&computations.V_norm,false);


        }

        if(gui.ground_truth.isChecked())
        {
            computations.GT=compute_ground_truth(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
            computations.GT_norm=computations.GT;
            computations.GT_norm.normalize();
            gui.canvas.push_obj(&computations.GT_norm,false);

        }
        if (gui.show_heat_map.isChecked())
        {

            if (gui.scalar_field.isChecked())
            {
                gui.canvas.pop(&computations.dual_m);
                gui.canvas.pop(&computations.m_grid);
                computations.dual_m.clear();


                computations.f=get_scalar_field(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());
                find_max_min_values(computations.f,computations.max,computations.min);
                computations.f_norm=heat_map_normalization(computations.f,computations.min,computations.max);
                computations.f_norm.copy_to_mesh(computations.m);
                computations.m.show_texture1D(TEXTURE_1D_HSV);


            }
            if (gui.Err.isChecked())
            {

                gui.canvas.pop(&computations.dual_m);
                gui.canvas.pop(&computations.m_grid);
                computations.dual_m.clear();

                if(gui.InsideError.isChecked())
                {

                    computations.f=get_scalar_field(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());
                    computations.V=compute_field(computations.m,computations.f,gui.method.currentIndex());
                    computations.GT_grid=compute_values_on_grid(computations.m_grid, gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());

                    bring_the_field_inside(computations.m,computations.m_grid,computations.V,computations.V_grid,gui.method.currentIndex());
                    computations.GT=compute_ground_truth(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
                    computations.err=estimate_error(computations.GT_grid,computations.V_grid,computations.m_grid,gui.mode.currentIndex(),LSDD,gui.ErrType.currentIndex(),gui.type_of_vertices.currentIndex());
                    computations.max=find_max_norm(computations.GT);
                    computations.err_norm=heat_map_normalization(computations.err,0,computations.max,gui.sl_error_neg.value(),gui.sl_error_pos.value());
                    computations.err_norm.copy_to_mesh(computations.m_grid);

                    gui.canvas.pop(&computations.m);
                    gui.canvas.push_obj(&computations.m_grid,false);
                    computations.m_grid.show_wireframe(false);
                    computations.m_grid.show_texture1D(TEXTURE_1D_HSV);
                    gui.canvas.push_obj(&computations.m,false);


                }else{


                    computations.f=get_scalar_field(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());
                    computations.GT=compute_ground_truth(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
                    computations.V=compute_field(computations.m,computations.f,gui.method.currentIndex());
                    computations.err=estimate_error(computations.GT,computations.V,computations.m,gui.mode.currentIndex(),gui.method.currentIndex(),gui.ErrType.currentIndex(),gui.type_of_vertices.currentIndex());
                    computations.max=find_max_norm(computations.GT);
                    computations.err_norm=heat_map_normalization(computations.err,0,computations.max,gui.sl_error_neg.value(),gui.sl_error_pos.value());

                    if(gui.method.currentIndex()==0)
                    {
                        computations.err_norm.copy_to_mesh(computations.dual_m);
                        gui.canvas.pop(&computations.m);
                        gui.canvas.push_obj(&computations.dual_m,false);
                        computations.dual_m.show_wireframe(false);
                        computations.dual_m.show_texture1D(TEXTURE_1D_HSV);
                        gui.canvas.push_obj(&computations.m,false);


                    }else
                    {
                        computations.err_norm.copy_to_mesh(computations.m);
                        computations.m.show_texture1D(TEXTURE_1D_HSV);

                    }
                }
            }
        }
        gui.canvas.updateGL();
    });


    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    QSpinBox::connect(&gui.a,static_cast<void(QSpinBox::*)(int)>(&QSpinBox::valueChanged), [&]()
    {
        computations.f=get_scalar_field(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());
        computations.GT=compute_ground_truth(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());

        if(gui.but.isChecked())
        {
            computations.f=get_scalar_field(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());
            computations.V=compute_field(computations.m,computations.f,gui.method.currentIndex());
            computations.V_norm=computations.V;
            computations.V_norm.normalize();
            gui.canvas.push_obj(&computations.V_norm,false);


        }
        if(gui.ground_truth.isChecked())
        {
            computations.GT=compute_ground_truth(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
            computations.GT_norm=computations.GT;
            computations.GT_norm.normalize();
            gui.canvas.push_obj(&computations.GT_norm,false);

        }

        if (gui.show_heat_map.isChecked())
        {
            if (gui.scalar_field.isChecked())
            {
                gui.canvas.pop(&computations.dual_m);
                gui.canvas.pop(&computations.m_grid);
                computations.dual_m.clear();


                computations.f=get_scalar_field(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());
                find_max_min_values(computations.f,computations.max,computations.min);
                computations.f_norm=heat_map_normalization(computations.f,computations.min,computations.max);
                computations.f_norm.copy_to_mesh(computations.m);
                computations.m.show_texture1D(TEXTURE_1D_HSV);


            }
            if (gui.Err.isChecked())
            {

                gui.canvas.pop(&computations.dual_m);
                gui.canvas.pop(&computations.m_grid);
                computations.dual_m.clear();

                if(gui.InsideError.isChecked())
                {

                    computations.f=get_scalar_field(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());
                    computations.V=compute_field(computations.m,computations.f,gui.method.currentIndex());
                    computations.GT_grid=compute_values_on_grid(computations.m_grid, gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());

                    bring_the_field_inside(computations.m,computations.m_grid,computations.V,computations.V_grid,gui.method.currentIndex());
                    computations.GT=compute_ground_truth(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
                    computations.err=estimate_error(computations.GT_grid,computations.V_grid,computations.m_grid,gui.mode.currentIndex(),LSDD,gui.ErrType.currentIndex(),gui.type_of_vertices.currentIndex());
                    computations.max=find_max_norm(computations.GT);
                    computations.err_norm=heat_map_normalization(computations.err,0,computations.max,gui.sl_error_neg.value(),gui.sl_error_pos.value());
                    computations.err_norm.copy_to_mesh(computations.m_grid);

                    gui.canvas.pop(&computations.m);
                    gui.canvas.push_obj(&computations.m_grid,false);
                    computations.m_grid.show_wireframe(false);
                    computations.m_grid.show_texture1D(TEXTURE_1D_HSV);
                    gui.canvas.push_obj(&computations.m,false);


                }else{


                    computations.f=get_scalar_field(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());
                    computations.GT=compute_ground_truth(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
                    computations.V=compute_field(computations.m,computations.f,gui.method.currentIndex());
                    computations.err=estimate_error(computations.GT,computations.V,computations.m,gui.mode.currentIndex(),gui.method.currentIndex(),gui.ErrType.currentIndex(),gui.type_of_vertices.currentIndex());
                    computations.max=find_max_norm(computations.GT);
                    computations.err_norm=heat_map_normalization(computations.err,0,computations.max,gui.sl_error_neg.value(),gui.sl_error_pos.value());

                    if(gui.method.currentIndex()==0)
                    {
                        computations.err_norm.copy_to_mesh(computations.dual_m);
                        gui.canvas.pop(&computations.m);
                        gui.canvas.push_obj(&computations.dual_m,false);
                        computations.dual_m.show_wireframe(false);
                        computations.dual_m.show_texture1D(TEXTURE_1D_HSV);
                        gui.canvas.push_obj(&computations.m,false);


                    }else
                    {
                        computations.err_norm.copy_to_mesh(computations.m);
                        computations.m.show_texture1D(TEXTURE_1D_HSV);

                    }
                }}

        }
        gui.canvas.updateGL();




    });
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    QSpinBox::connect(&gui.b,static_cast<void(QSpinBox::*)(int)>(&QSpinBox::valueChanged), [&]()
    {
        computations.f=get_scalar_field(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());
        computations.GT=compute_ground_truth(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());

        if(gui.but.isChecked())
        {
            computations.f=get_scalar_field(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());
            computations.V=compute_field(computations.m,computations.f,gui.method.currentIndex());
            computations.V_norm=computations.V;
            computations.V_norm.normalize();
            gui.canvas.push_obj(&computations.V_norm,false);


        }
        if(gui.ground_truth.isChecked())
        {
            computations.GT=compute_ground_truth(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
            computations.GT_norm=computations.GT;
            computations.GT_norm.normalize();
            gui.canvas.push_obj(&computations.GT_norm,false);

        }
        if (gui.show_heat_map.isChecked())
        {

            if (gui.scalar_field.isChecked())
            {
                gui.canvas.pop(&computations.dual_m);
                gui.canvas.pop(&computations.m_grid);
                computations.dual_m.clear();


                computations.f=get_scalar_field(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());
                find_max_min_values(computations.f,computations.max,computations.min);
                computations.f_norm=heat_map_normalization(computations.f,computations.min,computations.max);
                computations.f_norm.copy_to_mesh(computations.m);
                computations.m.show_texture1D(TEXTURE_1D_HSV);


            }
            if (gui.Err.isChecked())
            {

                gui.canvas.pop(&computations.dual_m);
                gui.canvas.pop(&computations.m_grid);
                computations.dual_m.clear();

                if(gui.InsideError.isChecked())
                {

                    computations.f=get_scalar_field(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());
                    computations.V=compute_field(computations.m,computations.f,gui.method.currentIndex());
                    computations.GT_grid=compute_values_on_grid(computations.m_grid, gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());

                    bring_the_field_inside(computations.m,computations.m_grid,computations.V,computations.V_grid,gui.method.currentIndex());
                    computations.GT=compute_ground_truth(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
                    computations.err=estimate_error(computations.GT_grid,computations.V_grid,computations.m_grid,gui.mode.currentIndex(),LSDD,gui.ErrType.currentIndex(),gui.type_of_vertices.currentIndex());
                    computations.max=find_max_norm(computations.GT);
                    computations.err_norm=heat_map_normalization(computations.err,0,computations.max,gui.sl_error_neg.value(),gui.sl_error_pos.value());
                    computations.err_norm.copy_to_mesh(computations.m_grid);

                    gui.canvas.pop(&computations.m);
                    gui.canvas.push_obj(&computations.m_grid,false);
                    computations.m_grid.show_wireframe(false);
                    computations.m_grid.show_texture1D(TEXTURE_1D_HSV);
                    gui.canvas.push_obj(&computations.m,false);


                }else{


                    computations.f=get_scalar_field(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());
                    computations.GT=compute_ground_truth(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
                    computations.V=compute_field(computations.m,computations.f,gui.method.currentIndex());
                    computations.err=estimate_error(computations.GT,computations.V,computations.m,gui.mode.currentIndex(),gui.method.currentIndex(),gui.ErrType.currentIndex(),gui.type_of_vertices.currentIndex());
                    computations.max=find_max_norm(computations.GT);
                    computations.err_norm=heat_map_normalization(computations.err,0,computations.max,gui.sl_error_neg.value(),gui.sl_error_pos.value());

                    if(gui.method.currentIndex()==0)
                    {
                        computations.err_norm.copy_to_mesh(computations.dual_m);
                        gui.canvas.pop(&computations.m);
                        gui.canvas.push_obj(&computations.dual_m,false);
                        computations.dual_m.show_wireframe(false);
                        computations.dual_m.show_texture1D(TEXTURE_1D_HSV);
                        gui.canvas.push_obj(&computations.m,false);


                    }else
                    {
                        computations.err_norm.copy_to_mesh(computations.m);
                        computations.m.show_texture1D(TEXTURE_1D_HSV);

                    }
                }}
        }

        gui.canvas.updateGL();
    });
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    QSpinBox::connect(&gui.c,static_cast<void(QSpinBox::*)(int)>(&QSpinBox::valueChanged), [&]()
    {
        computations.f=get_scalar_field(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());
        computations.GT=compute_ground_truth(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());

        if(gui.but.isChecked())
        {
            computations.f=get_scalar_field(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());
            computations.V=compute_field(computations.m,computations.f,gui.method.currentIndex());
            computations.V_norm=computations.V;
            computations.V_norm.normalize();
            gui.canvas.push_obj(&computations.V_norm,false);


        }
        if(gui.ground_truth.isChecked())
        {
            computations.GT=compute_ground_truth(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
            computations.GT_norm=computations.GT;
            computations.GT_norm.normalize();
            gui.canvas.push_obj(&computations.GT_norm,false);
            gui.canvas.updateGL();
        }
        if (gui.show_heat_map.isChecked())
        {
            if (gui.scalar_field.isChecked())
            {
                gui.canvas.pop(&computations.dual_m);
                gui.canvas.pop(&computations.m_grid);
                computations.dual_m.clear();


                computations.f=get_scalar_field(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());
                find_max_min_values(computations.f,computations.max,computations.min);
                computations.f_norm=heat_map_normalization(computations.f,computations.min,computations.max);
                computations.f_norm.copy_to_mesh(computations.m);
                computations.m.show_texture1D(TEXTURE_1D_HSV);


            }
            if (gui.Err.isChecked())
            {

                gui.canvas.pop(&computations.dual_m);
                gui.canvas.pop(&computations.m_grid);
                computations.dual_m.clear();

                if(gui.InsideError.isChecked())
                {

                    computations.f=get_scalar_field(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());
                    computations.V=compute_field(computations.m,computations.f,gui.method.currentIndex());
                    computations.GT_grid=compute_values_on_grid(computations.m_grid, gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());

                    bring_the_field_inside(computations.m,computations.m_grid,computations.V,computations.V_grid,gui.method.currentIndex());
                    computations.GT=compute_ground_truth(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
                    computations.err=estimate_error(computations.GT_grid,computations.V_grid,computations.m_grid,gui.mode.currentIndex(),LSDD,gui.ErrType.currentIndex(),gui.type_of_vertices.currentIndex());
                    computations.max=find_max_norm(computations.GT);
                    computations.err_norm=heat_map_normalization(computations.err,0,computations.max,gui.sl_error_neg.value(),gui.sl_error_pos.value());
                    computations.err_norm.copy_to_mesh(computations.m_grid);

                    gui.canvas.pop(&computations.m);
                    gui.canvas.push_obj(&computations.m_grid,false);
                    computations.m_grid.show_wireframe(false);
                    computations.m_grid.show_texture1D(TEXTURE_1D_HSV);
                    gui.canvas.push_obj(&computations.m,false);

                }else{


                    computations.f=get_scalar_field(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());
                    computations.GT=compute_ground_truth(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
                    computations.V=compute_field(computations.m,computations.f,gui.method.currentIndex());
                    computations.err=estimate_error(computations.GT,computations.V,computations.m,gui.mode.currentIndex(),gui.method.currentIndex(),gui.ErrType.currentIndex(),gui.type_of_vertices.currentIndex());
                    computations.max=find_max_norm(computations.GT);
                    computations.err_norm=heat_map_normalization(computations.err,0,computations.max,gui.sl_error_neg.value(),gui.sl_error_pos.value());

                    if(gui.method.currentIndex()==0)
                    {
                        computations.err_norm.copy_to_mesh(computations.dual_m);
                        gui.canvas.pop(&computations.m);
                        gui.canvas.push_obj(&computations.dual_m,false);
                        computations.dual_m.show_wireframe(false);
                        computations.dual_m.show_texture1D(TEXTURE_1D_HSV);
                        gui.canvas.push_obj(&computations.m,false);


                    }else
                    {
                        computations.err_norm.copy_to_mesh(computations.m);
                        computations.m.show_texture1D(TEXTURE_1D_HSV);

                    }
                }}
        }

        gui.canvas.updateGL();

    });
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    QRadioButton::connect(&gui.scalar_field, &QRadioButton::toggled, [&]()
    {


        if (gui.scalar_field.isChecked())
        {
            gui.canvas.pop(&computations.dual_m);
            gui.canvas.pop(&computations.m_grid);

            gui.error.setDisabled(true);
            computations.f=get_scalar_field(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());
            find_max_min_values(computations.f,computations.max,computations.min);
            computations.f_norm=heat_map_normalization(computations.f,computations.min,computations.max);

            computations.f_norm.copy_to_mesh(computations.m);
            computations.m.show_texture1D(TEXTURE_1D_HSV);

            if(gui.but.isChecked())
            {
                computations.f=get_scalar_field(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());
                computations.V=compute_field(computations.m,computations.f,gui.method.currentIndex());
                computations.V_norm=computations.V;
                computations.V_norm.normalize();
                gui.canvas.push_obj(&computations.V_norm,false);


            }
            if(gui.ground_truth.isChecked())
            {
                computations.GT=compute_ground_truth(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
                computations.GT_norm=computations.GT;
                computations.GT_norm.normalize();
                gui.canvas.push_obj(&computations.GT_norm,false);

            }




        }
        gui.canvas.updateGL();


    });
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    QRadioButton::connect(&gui.Err, &QRadioButton::toggled, [&]()
    {

        if (gui.Err.isChecked())
        {
            gui.error.setDisabled(false);
            gui.canvas.pop(&computations.dual_m);
            gui.canvas.pop(&computations.m_grid);


            if(gui.InsideError.isChecked())
            {
                computations.f=get_scalar_field(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());
                computations.V=compute_field(computations.m,computations.f,gui.method.currentIndex());
                computations.GT_grid=compute_values_on_grid(computations.m_grid, gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());

                bring_the_field_inside(computations.m,computations.m_grid,computations.V,computations.V_grid,gui.method.currentIndex());
                computations.GT=compute_ground_truth(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
                computations.err=estimate_error(computations.GT_grid,computations.V_grid,computations.m_grid,gui.mode.currentIndex(),LSDD,gui.ErrType.currentIndex(),gui.type_of_vertices.currentIndex());
                computations.max=find_max_norm(computations.GT);
                computations.err_norm=heat_map_normalization(computations.err,0,computations.max,gui.sl_error_neg.value(),gui.sl_error_pos.value());
                computations.err_norm.copy_to_mesh(computations.m_grid);

                gui.canvas.pop(&computations.m);
                gui.canvas.push_obj(&computations.m_grid,false);
                computations.m_grid.show_wireframe(false);
                computations.m_grid.show_texture1D(TEXTURE_1D_HSV);
                gui.canvas.push_obj(&computations.m,false);

            }else{




                computations.f=get_scalar_field(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());
                computations.GT=compute_ground_truth(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
                computations.V=compute_field(computations.m,computations.f,gui.method.currentIndex());
                computations.err=estimate_error(computations.GT,computations.V,computations.m,gui.mode.currentIndex(),gui.method.currentIndex(),gui.ErrType.currentIndex(),gui.type_of_vertices.currentIndex());
                computations.max=find_max_norm(computations.GT);
                computations.err_norm=heat_map_normalization(computations.err,0,computations.max,gui.sl_error_neg.value(),gui.sl_error_pos.value());

                if(gui.method.currentIndex()==0)
                {
                    computations.err_norm.copy_to_mesh(computations.dual_m);
                    gui.canvas.pop(&computations.m);
                    gui.canvas.push_obj(&computations.dual_m,false);
                    computations.dual_m.show_wireframe(false);
                    computations.dual_m.show_texture1D(TEXTURE_1D_HSV);
                    gui.canvas.push_obj(&computations.m,false);



                }else
                {
                    computations.err_norm.copy_to_mesh(computations.m);
                    computations.m.show_texture1D(TEXTURE_1D_HSV);


                }
            }
            if(gui.but.isChecked())
            {
                computations.f=get_scalar_field(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());
                computations.V=compute_field(computations.m,computations.f,gui.method.currentIndex());
                computations.V_norm=computations.V;
                computations.V_norm.normalize();
                gui.canvas.push_obj(&computations.V_norm,false);


            }
            if(gui.ground_truth.isChecked())
            {
                computations.GT=compute_ground_truth(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
                computations.GT_norm=computations.GT;
                computations.GT_norm.normalize();
                gui.canvas.push_obj(&computations.GT_norm,false);

            }

        }
        gui.canvas.updateGL();


    });
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    QCheckBox::connect(&gui.InsideError, &QCheckBox::stateChanged, [&]()
    {
        gui.canvas.pop(&computations.dual_m);
        gui.canvas.pop(&computations.m_grid);
        computations.dual_m.clear();


        if(gui.InsideError.isChecked())
        {


                computations.f=get_scalar_field(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());
                computations.V=compute_field(computations.m,computations.f,gui.method.currentIndex());
                computations.GT_grid=compute_values_on_grid(computations.m_grid, gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());

                bring_the_field_inside(computations.m,computations.m_grid,computations.V,computations.V_grid,gui.method.currentIndex());
                computations.GT=compute_ground_truth(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
                computations.err=estimate_error(computations.GT_grid,computations.V_grid,computations.m_grid,gui.mode.currentIndex(),LSDD,gui.ErrType.currentIndex(),gui.type_of_vertices.currentIndex());
                computations.max=find_max_norm(computations.GT);
                computations.err_norm=heat_map_normalization(computations.err,0,computations.max,gui.sl_error_neg.value(),gui.sl_error_pos.value());
                computations.err_norm.copy_to_mesh(computations.m_grid);

                gui.canvas.pop(&computations.m);
                gui.canvas.push_obj(&computations.m_grid,false);
                computations.m_grid.show_wireframe(false);
                computations.m_grid.show_texture1D(TEXTURE_1D_HSV);
                gui.canvas.push_obj(&computations.m,false);

                if(gui.but.isChecked())
                {
                    computations.f=get_scalar_field(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());
                    computations.V=compute_field(computations.m,computations.f,gui.method.currentIndex());
                    computations.V_norm=computations.V;
                    computations.V_norm.normalize();
                    gui.canvas.push_obj(&computations.V_norm,false);


                }
                if(gui.ground_truth.isChecked())
                {
                    computations.GT=compute_ground_truth(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
                    computations.GT_norm=computations.GT;
                    computations.GT_norm.normalize();
                    gui.canvas.push_obj(&computations.GT_norm,false);

                }

        }else{


            computations.f=get_scalar_field(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());
            computations.GT=compute_ground_truth(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
            computations.V=compute_field(computations.m,computations.f,gui.method.currentIndex());
            computations.err=estimate_error(computations.GT,computations.V,computations.m,gui.mode.currentIndex(),gui.method.currentIndex(),gui.ErrType.currentIndex(),gui.type_of_vertices.currentIndex());
            computations.max=find_max_norm(computations.GT);
            computations.err_norm=heat_map_normalization(computations.err,0,computations.max,gui.sl_error_neg.value(),gui.sl_error_pos.value());

            if(gui.method.currentIndex()==0)
            {
                computations.err_norm.copy_to_mesh(computations.dual_m);
                gui.canvas.pop(&computations.m);
                gui.canvas.push_obj(&computations.dual_m,false);
                computations.dual_m.show_wireframe(false);
                computations.dual_m.show_texture1D(TEXTURE_1D_HSV);
                gui.canvas.push_obj(&computations.m,false);



            }else
            {
                computations.err_norm.copy_to_mesh(computations.m);
                computations.m.show_texture1D(TEXTURE_1D_HSV);


            }





        }

        gui.canvas.updateGL();
    });
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    QComboBox::connect(&gui.type_of_vertices, static_cast<void(QComboBox::*)(int)>(&QComboBox::activated), [&]()
    {
        gui.canvas.pop(&computations.dual_m);
        gui.canvas.pop(&computations.m_grid);
        if(gui.InsideError.isChecked())
        {

            computations.f=get_scalar_field(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());
            computations.GT_grid=compute_values_on_grid(computations.m_grid, gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());
            bring_the_field_inside(computations.m,computations.m_grid,computations.V,computations.V_grid,gui.method.currentIndex());
            computations.GT=compute_ground_truth(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
            computations.err=estimate_error(computations.GT_grid,computations.V_grid,computations.m_grid,gui.mode.currentIndex(),LSDD,gui.ErrType.currentIndex(),gui.type_of_vertices.currentIndex());
            computations.max=find_max_norm(computations.GT);
            computations.err_norm=heat_map_normalization(computations.err,0,computations.max,gui.sl_error_neg.value(),gui.sl_error_pos.value());
            computations.err_norm.copy_to_mesh(computations.m_grid);

            gui.canvas.pop(&computations.m);
            gui.canvas.push_obj(&computations.m_grid,false);
            computations.m_grid.show_wireframe(false);
            computations.m_grid.show_texture1D(TEXTURE_1D_HSV);
            gui.canvas.push_obj(&computations.m,false);



        }else

        {
            computations.f=get_scalar_field(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());
            computations.GT=compute_ground_truth(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
            computations.V=compute_field(computations.m,computations.f,gui.method.currentIndex());
            computations.err=estimate_error(computations.GT,computations.V,computations.m,gui.mode.currentIndex(),gui.method.currentIndex(),gui.ErrType.currentIndex(),gui.type_of_vertices.currentIndex());
            computations.max=find_max_norm(computations.GT);
            computations.err_norm=heat_map_normalization(computations.err,0,computations.max,gui.sl_error_neg.value(),gui.sl_error_pos.value());

            if(gui.method.currentIndex()==0)
            {
                computations.err_norm.copy_to_mesh(computations.dual_m);
                gui.canvas.pop(&computations.m);
                gui.canvas.push_obj(&computations.dual_m,false);
                computations.dual_m.show_wireframe(false);
                computations.dual_m.show_texture1D(TEXTURE_1D_HSV);
                gui.canvas.push_obj(&computations.m,false);


            }else
            {
                computations.err_norm.copy_to_mesh(computations.m);
                computations.m.show_texture1D(TEXTURE_1D_HSV);


            }
        }
        gui.canvas.updateGL();

    });
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    QComboBox::connect(&gui.ErrType, static_cast<void(QComboBox::*)(int)>(&QComboBox::activated), [&]()
    {

        gui.canvas.pop(&computations.dual_m);
        gui.canvas.pop(&computations.m_grid);
        computations.dual_m.clear();


        if(gui.InsideError.isChecked())
        {

            computations.f=get_scalar_field(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());
            computations.GT_grid=compute_values_on_grid(computations.m_grid, gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());
            bring_the_field_inside(computations.m,computations.m_grid,computations.V,computations.V_grid,gui.method.currentIndex());
            computations.GT=compute_ground_truth(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
            computations.err=estimate_error(computations.GT_grid,computations.V_grid,computations.m_grid,gui.mode.currentIndex(),LSDD,gui.ErrType.currentIndex(),gui.type_of_vertices.currentIndex());
            computations.max=find_max_norm(computations.GT);
            computations.err_norm=heat_map_normalization(computations.err,0,computations.max,gui.sl_error_neg.value(),gui.sl_error_pos.value());
            computations.err_norm.copy_to_mesh(computations.m_grid);

            gui.canvas.pop(&computations.m);
            gui.canvas.push_obj(&computations.m_grid,false);
            computations.m_grid.show_wireframe(false);
            computations.m_grid.show_texture1D(TEXTURE_1D_HSV);
            gui.canvas.push_obj(&computations.m,false);



        }else{

            computations.f=get_scalar_field(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());
            computations.GT=compute_ground_truth(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
            computations.V=compute_field(computations.m,computations.f,gui.method.currentIndex());
            computations.err=estimate_error(computations.GT,computations.V,computations.m,gui.mode.currentIndex(),gui.method.currentIndex(),gui.ErrType.currentIndex(),gui.type_of_vertices.currentIndex());
            computations.max=find_max_norm(computations.GT);
            computations.err_norm=heat_map_normalization(computations.err,0,computations.max,gui.sl_error_neg.value(),gui.sl_error_pos.value());

            if(gui.method.currentIndex()==0)
            {
                computations.err_norm.copy_to_mesh(computations.dual_m);
                gui.canvas.pop(&computations.m);
                gui.canvas.push_obj(&computations.dual_m,false);
                computations.dual_m.show_wireframe(false);
                computations.dual_m.show_texture1D(TEXTURE_1D_HSV);
                gui.canvas.push_obj(&computations.m,false);


            }else
            {
                computations.err_norm.copy_to_mesh(computations.m);
                computations.m.show_texture1D(TEXTURE_1D_HSV);


            }

        }

        if(gui.but.isChecked())
        {
            computations.f=get_scalar_field(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());
            computations.V=compute_field(computations.m,computations.f,gui.method.currentIndex());
            computations.V_norm=computations.V;
            computations.V_norm.normalize();
            gui.canvas.push_obj(&computations.V_norm,false);


        }
        if(gui.ground_truth.isChecked())
        {
            computations.GT=compute_ground_truth(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
            computations.GT_norm=computations.GT;
            computations.GT_norm.normalize();
            gui.canvas.push_obj(&computations.GT_norm,false);

        }


        gui.canvas.updateGL();

    });
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    QComboBox::connect(&gui.mode, static_cast<void(QComboBox::*)(int)>(&QComboBox::activated), [&]()
    {


        gui.canvas.pop(&computations.dual_m);
        gui.canvas.pop(&computations.m_grid);
        computations.dual_m.clear();


        if(gui.InsideError.isChecked())
        {
            computations.f=get_scalar_field(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());
            computations.GT_grid=compute_values_on_grid(computations.m_grid, gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());
            bring_the_field_inside(computations.m,computations.m_grid,computations.V,computations.V_grid,gui.method.currentIndex());
            computations.GT=compute_ground_truth(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
            computations.err=estimate_error(computations.GT_grid,computations.V_grid,computations.m_grid,gui.mode.currentIndex(),LSDD,gui.ErrType.currentIndex(),gui.type_of_vertices.currentIndex());
            computations.max=find_max_norm(computations.GT);
            computations.err_norm=heat_map_normalization(computations.err,0,computations.max,gui.sl_error_neg.value(),gui.sl_error_pos.value());
            computations.err_norm.copy_to_mesh(computations.m_grid);

            gui.canvas.pop(&computations.m);
            gui.canvas.push_obj(&computations.m_grid,false);
            computations.m_grid.show_wireframe(false);
            computations.m_grid.show_texture1D(TEXTURE_1D_HSV);
            gui.canvas.push_obj(&computations.m,false);


        }else{


            computations.f=get_scalar_field(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());
            computations.GT=compute_ground_truth(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
            computations.V=compute_field(computations.m,computations.f,gui.method.currentIndex());
            computations.err=estimate_error(computations.GT,computations.V,computations.m,gui.mode.currentIndex(),gui.method.currentIndex(),gui.ErrType.currentIndex(),gui.type_of_vertices.currentIndex());
            computations.max=find_max_norm(computations.GT);
            computations.err_norm=heat_map_normalization(computations.err,0,computations.max,gui.sl_error_neg.value(),gui.sl_error_pos.value());

            if(gui.method.currentIndex()==0)
            {
                computations.err_norm.copy_to_mesh(computations.dual_m);
                gui.canvas.pop(&computations.m);
                gui.canvas.push_obj(&computations.dual_m,false);
                computations.dual_m.show_wireframe(false);
                computations.dual_m.show_texture1D(TEXTURE_1D_HSV);
                gui.canvas.push_obj(&computations.m,false);


            }else
            {
                computations.err_norm.copy_to_mesh(computations.m);
                computations.m.show_texture1D(TEXTURE_1D_HSV);


            }
        }

        if(gui.but.isChecked())
        {
            computations.f=get_scalar_field(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());
            computations.V=compute_field(computations.m,computations.f,gui.method.currentIndex());
            computations.V_norm=computations.V;
            computations.V_norm.normalize();
            gui.canvas.push_obj(&computations.V_norm,false);


        }
        if(gui.ground_truth.isChecked())
        {
            computations.GT=compute_ground_truth(computations.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
            computations.GT_norm=computations.GT;
            computations.GT_norm.normalize();
            gui.canvas.push_obj(&computations.GT_norm,false);

        }


        gui.canvas.updateGL();
    });

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    QPushButton::connect(&gui.save, &QPushButton::clicked, [&]()
    {
//        ScalarField data=computations.err;
//        std::ofstream outfile;


//        outfile.open("/Users/Dala/Documents/MATLAB/Gradient_fields/" + std::to_string(10) + "result.csv");
//        for(int i=0; i<computations.m.num_verts();++i)
//        {
//            outfile<<data[i]<<","<<std::endl;
//        }


                std::string filename = QFileDialog::getSaveFileName(NULL, "Save mesh", ".", "3D Meshes (*.off *.obj *.iv);; OBJ(*.obj);; OFF(*.off);; IV(*.iv)").toStdString();
                if (!filename.empty()) computations.m.save(filename.c_str());




    });
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    QPushButton::connect(&gui.reset, &QPushButton::clicked, [&]()
    {

        gui.show_heat_map.setChecked(false);
        gui.choose_heat_map.setDisabled(true);

        gui.sl_error_neg.setDisabled(true);
        gui.sl_error_pos.setDisabled(true);

        gui.error.setDisabled(true);

        gui.but.setChecked(false);


        gui.ground_truth.setChecked(false);

        computations.m.show_vert_color();
        gui.canvas.pop_all_occurrences_of(DRAWABLE_VECTOR_FIELD);
        gui.canvas.pop(&computations.dual_m);
        gui.canvas.updateGL();

    });
}






