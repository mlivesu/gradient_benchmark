#include <QGridLayout>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QLabel>
#include "gui.h"
#include"functions_cubic.h"
#include <QDir>
#include<iostream>
#include<fstream>



extern COMPUTATIONS computations_cubic;

using namespace cinolib;

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void init_gui(GUI & gui)
{
    QGridLayout *global_layout  = new QGridLayout();
    QGridLayout *choose_heat_map_layout=new QGridLayout();
    QGroupBox *show= new QGroupBox;
    QGroupBox *slicing=new QGroupBox;
    QGridLayout *error_layout=new QGridLayout();
    QGridLayout *slicing_layout=new QGridLayout();
    QGridLayout *show_layout=new QGridLayout();



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
    button_layout->addWidget(new QLabel("Tetrahedralization:"),3,0);
    button_layout->addWidget(&gui.choose_tri,3,1);
    gui.choose_tri.setMaximumSize(150,50);
    button_layout->addWidget(new QLabel("Number of samples:"),4,0);
    button_layout->addWidget(&gui.N,4,1);
    gui.N.setMaximumSize(100,50);
    gui.N.setMinimum(1);
    gui.N.setMaximum(500);
    gui.N.setValue(10);

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
    gui.error.setTitle("Error");
    gui.error.setLayout(error_layout);
    gui.error.setDisabled(true);

    button_layout->addWidget(&gui.error,12,1);

    slicing_layout->addWidget(new QLabel("X"),0,0);
    slicing_layout->addWidget(&gui.sl_slice_x,0,1);
    gui.sl_slice_x.setOrientation(Qt::Horizontal);
    slicing_layout->addWidget(new QLabel("Flip"),0,2);
    slicing_layout->addWidget(&gui.cb_slice_flip_x,0,3);

    slicing_layout->addWidget(new QLabel("Y"),1,0);
    slicing_layout->addWidget(&gui.sl_slice_y,1,1);
    gui.sl_slice_y.setOrientation(Qt::Horizontal);
    slicing_layout->addWidget(new QLabel("Flip"),1,2);
    slicing_layout->addWidget(&gui.cb_slice_flip_y,1,3);

    slicing_layout->addWidget(new QLabel("Z"),2,0);
    slicing_layout->addWidget(&gui.sl_slice_z,2,1);
    gui.sl_slice_z.setOrientation(Qt::Horizontal);
    slicing_layout->addWidget(new QLabel("Flip"),2,2);
    slicing_layout->addWidget(&gui.cb_slice_flip_z,2,3);

    gui.sl_slice_x.setValue(gui.sl_slice_x.maximum());
    gui.sl_slice_y.setValue(gui.sl_slice_y.maximum());
    gui.sl_slice_z.setValue(gui.sl_slice_z.maximum());

    gui.slicing.setTitle("Slicing");
    gui.slicing.setLayout(slicing_layout);

    button_layout->addWidget(&gui.slicing,13,1);





    button_layout->addWidget(&gui.save,16,0);
    button_layout->addWidget(&gui.reset,16,1);
     button_layout->addWidget(&gui.load,16,2);
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


    button_layout->setSpacing(0);
    button_layout->setHorizontalSpacing(0);
    button_layout->setContentsMargins(0,0,0,0);
    toolbar_layout->addLayout(button_layout);
    toolbar_layout->setSpacing(0);

    toolbar_layout->setContentsMargins(0,0,0,0);


    //

    for(int i=0; i<N_FUNCTIONS; ++i) gui.cb_function.addItem(f_names[i].c_str());
    for(int i=0; i<N_METHODS; ++i) gui.method.addItem(method_names[i].c_str());
    for(int i=0; i<N_tri; ++i) gui.choose_tri.addItem(tri_names[i].c_str());
    for(int i=0; i<=relative; ++i) gui.ErrType.addItem(err_names[i].c_str());
    for(int i=0; i<=angle; ++i) gui.mode.addItem(mode_names[i].c_str());
    for(int i=0; i<=only_srf; ++i) gui.type_of_vertices.addItem(vertices_names[i].c_str());
    gui.but.setText("Gradient");
    gui.ground_truth.setText("Ground Truth");
    gui.wireframe.setText("Wireframe");
    gui.scalar_field.setText("Scalar Field");
    gui.Err.setText("Error");
    gui.save.setText("Save");
    gui.load.setText("Load");
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
                uint vid;
                if(gui.InsideError.isChecked())
                {
                    vid=closest_vertex(p,&computations_cubic.m_grid);
                }else if (gui.Err.isChecked())
                {
                    if(gui.method.currentIndex()!=0)
                    {
                        vid= closest_vertex(p,&computations_cubic.m);
                    }else
                    {
                        vid = closest_vertex(p,&computations_cubic.dual_m);
                    }

                }
                std::cout << mode_names[gui.mode.currentIndex()].c_str() <<" error calculated:"<<std::endl;
                std::cout<<computations_cubic.err[vid] << std::endl;
            }
        }
        else if(e->modifiers() == Qt::ShiftModifier)
        {

            vec3d p;
            vec2i click(e->x(), e->y());
            if (c->unproject(click, p)) // transform click in a 3d point
            {

                uint vid = closest_vertex(p,&computations_cubic.m);
                std::cout <<"Scalar Field value:"<< computations.f[vid] << std::endl;

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

            gui.sl_error_neg.setDisabled(false);
            gui.sl_error_pos.setDisabled(false);

            if (gui.scalar_field.isChecked())
            {
                gui.canvas.pop(&computations_cubic.dual_m);
                gui.canvas.pop(&computations_cubic.m_grid);

                computations_cubic.f=get_scalar_field(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),0,0);
                find_max_min_values(computations_cubic.f,computations_cubic.max,computations_cubic.min);
                computations_cubic.f_norm=heat_map_normalization(computations_cubic.f,computations_cubic.min,computations_cubic.max);
                computations_cubic.f_norm.copy_to_mesh(computations_cubic.m);
                computations_cubic.m.show_out_texture1D(TEXTURE_1D_HSV);


            }

            if (gui.Err.isChecked())
            {
                gui.error.setDisabled(false);

                gui.canvas.pop(&computations_cubic.dual_m);
                gui.canvas.pop(&computations_cubic.m_grid);

                if(gui.InsideError.isChecked())
                {/*
                    computations_cubic.f=get_scalar_field(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),0,0);
                    computations_cubic.V=compute_field(computations_cubic.m,computations_cubic.f,gui.method.currentIndex());
                    //computations_cubic.GT_grid=compute_values_on_grid(computations_cubic.m_grid, gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());

                    bring_the_field_inside(computations_cubic.m,computations_cubic.m_grid,computations_cubic.V,computations_cubic.V_grid,gui.method.currentIndex());
                    computations_cubic.GT=compute_ground_truth(computations_cubic.m_grid,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),LSDD);
                    computations_cubic.err=estimate_error(computations_cubic.GT_grid,computations_cubic.V_grid,computations_cubic.m_grid,gui.mode.currentIndex(),LSDD,gui.ErrType.currentIndex(),gui.type_of_vertices.currentIndex());
                    computations_cubic.max=find_max_norm(computations_cubic.GT);*/
                    computations_cubic.err_norm=heat_map_normalization(computations_cubic.err,0,computations_cubic.max,gui.sl_error_neg.value(),gui.sl_error_pos.value());
                    computations_cubic.err_norm.copy_to_mesh(computations_cubic.m_grid);

                    gui.canvas.pop(&computations_cubic.m);
                    gui.canvas.push_obj(&computations_cubic.m_grid,false);
                    computations_cubic.m_grid.show_in_wireframe(false);
                    computations_cubic.m_grid.show_in_texture1D(TEXTURE_1D_HSV);
                    computations_cubic.m.show_in_vert_color();
                    computations_cubic.m.show_out_vert_color();
                    gui.canvas.push_obj(&computations_cubic.m,false);
                }else
                {
                    //                    computations_cubic.f=get_scalar_field(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),0,0);
                    //                    computations_cubic.GT=compute_ground_truth(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
                    //                    computations_cubic.V=compute_field(computations_cubic.m,computations_cubic.f,gui.method.currentIndex());
                    //                    computations_cubic.err=estimate_error(computations_cubic.GT_grid,computations_cubic.V_grid,computations_cubic.m_grid,gui.mode.currentIndex(),LSDD,gui.ErrType.currentIndex(),gui.type_of_vertices.currentIndex());
                    //                    computations_cubic.max=find_max_norm(computations_cubic.GT);
                    computations_cubic.err_norm=heat_map_normalization(computations_cubic.err,0,computations_cubic.max,gui.sl_error_neg.value(),gui.sl_error_pos.value());

                    if(gui.method.currentIndex()==0)
                    {


                        computations_cubic.err_norm.copy_to_mesh(computations_cubic.dual_m);

                        gui.canvas.pop(&computations_cubic.m);
                        gui.canvas.push_obj(&computations_cubic.dual_m,false);
                        computations_cubic.dual_m.show_in_wireframe(false);
                        computations_cubic.dual_m.show_out_wireframe(false);

                        computations_cubic.dual_m.show_in_texture1D(TEXTURE_1D_HSV);
                        gui.canvas.push_obj(&computations_cubic.m,false);
                    }else
                    {
                        computations_cubic.err_norm.copy_to_mesh(computations_cubic.m);
                        computations_cubic.m.show_in_texture1D(TEXTURE_1D_HSV);


                    }
                }



                if(gui.but.isChecked())
                {
                    //                    computations_cubic.f=get_scalar_field(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),0,0);
                    //                    computations_cubic.V=compute_field(computations_cubic.m,computations_cubic.f,gui.method.currentIndex());
                    //                    computations_cubic.V_norm=computations_cubic.V;
                    //                    computations_cubic.V_norm.normalize();
                    gui.canvas.push_obj(&computations_cubic.V_norm,false);


                }
                if(gui.ground_truth.isChecked())
                {
                    //                    computations_cubic.GT=compute_ground_truth(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
                    //                    computations_cubic.GT_norm=computations_cubic.GT;
                    //                    computations_cubic.GT_norm.normalize();
                    gui.canvas.push_obj(&computations_cubic.GT_norm,false);

                }

            }

        }else{

            gui.choose_heat_map.setDisabled(true);
            gui.sl_error_neg.setDisabled(true);
            gui.sl_error_pos.setDisabled(true);
            gui.error.setDisabled(true);
            computations_cubic.dual_m.show_in_vert_color();
            computations_cubic.dual_m_grid.show_in_vert_color();
            computations_cubic.dual_m.show_in_vert_color();
        }


        gui.canvas.updateGL();
    });
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    QCheckBox::connect(&gui.but, &QCheckBox::stateChanged, [&]()
    {
        if(gui.but.isChecked())
        {
                        computations_cubic.f=get_scalar_field(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),0,0);
                        computations_cubic.V=compute_field(computations_cubic.m,computations_cubic.f,gui.method.currentIndex());
                        computations_cubic.V_norm=computations_cubic.V;
                        computations_cubic.V_norm.normalize();
            gui.canvas.push_obj(&computations_cubic.V_norm,false);


        }
        else
        {
            gui.canvas.pop_all_occurrences_of(DRAWABLE_VECTOR_FIELD);

            if(gui.ground_truth.isChecked())
            {
                //                computations_cubic.GT=compute_ground_truth(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
                //                computations_cubic.GT_norm=computations_cubic.GT;
                //                computations_cubic.GT_norm.normalize();
                gui.canvas.push_obj(&computations_cubic.GT_norm,false);

            }
        }
        gui.canvas.updateGL();



    });
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    QSlider::connect(&gui.sl_error_neg, &QSlider::valueChanged, [&]()
    {


        if (gui.Err.isChecked())
        {
            gui.canvas.pop(&computations_cubic.dual_m);
            gui.canvas.pop(&computations_cubic.m_grid);
            if(gui.InsideError.isChecked())
            {
                computations_cubic.f=get_scalar_field(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),0,1);
                computations_cubic.V=compute_field(computations_cubic.m,computations_cubic.f,gui.method.currentIndex());


                bring_the_field_inside(computations_cubic.m,computations_cubic.m_grid,computations_cubic.V,computations_cubic.V_grid,gui.method.currentIndex());
                computations_cubic.GT=compute_ground_truth(computations_cubic.m_grid,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());

                computations_cubic.err=estimate_error(computations_cubic.GT,computations_cubic.V_grid,computations_cubic.m_grid,gui.mode.currentIndex(),LSDD,gui.ErrType.currentIndex(),gui.type_of_vertices.currentIndex());

                computations_cubic.max=find_max_norm(computations_cubic.GT);
                computations_cubic.err_norm=heat_map_normalization(computations_cubic.err,0,computations_cubic.max,gui.sl_error_neg.value(),gui.sl_error_pos.value());
                computations_cubic.err_norm.copy_to_mesh(computations_cubic.m_grid);

                gui.canvas.pop(&computations_cubic.m);
                gui.canvas.push_obj(&computations_cubic.m_grid,false);
                computations_cubic.m_grid.show_in_wireframe(false);
                computations_cubic.m_grid.show_out_wireframe(false);
                computations_cubic.m_grid.show_in_texture1D(TEXTURE_1D_HSV);
                computations_cubic.m_grid.show_out_texture1D(TEXTURE_1D_HSV);
                gui.canvas.push_obj(&computations_cubic.m,false);

            }else{

                computations_cubic.f=get_scalar_field(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),0,1);
                computations_cubic.GT=compute_ground_truth(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
                computations_cubic.V=compute_field(computations_cubic.m,computations_cubic.f,gui.method.currentIndex());
                computations_cubic.err=estimate_error(computations_cubic.GT,computations_cubic.V,computations_cubic.m,gui.mode.currentIndex(),gui.method.currentIndex(),gui.ErrType.currentIndex(),gui.type_of_vertices.currentIndex());
                computations_cubic.max=find_max_norm(computations_cubic.GT);
                computations_cubic.err_norm=heat_map_normalization(computations_cubic.err,0,computations_cubic.max,gui.sl_error_neg.value(),gui.sl_error_pos.value());

                if(gui.method.currentIndex()==0)
                {


                    computations_cubic.err_norm.copy_to_mesh(computations_cubic.dual_m);

                    gui.canvas.pop(&computations_cubic.m);
                    gui.canvas.push_obj(&computations_cubic.dual_m,false);
                    computations_cubic.dual_m.show_in_wireframe(false);
                    computations_cubic.dual_m.show_out_wireframe(false);

                    computations_cubic.dual_m.show_in_texture1D(TEXTURE_1D_HSV);
                    computations_cubic.dual_m.show_out_texture1D(TEXTURE_1D_HSV);
                    gui.canvas.push_obj(&computations_cubic.m,false);
                }else
                {
                    computations_cubic.err_norm.copy_to_mesh(computations_cubic.m);
                    computations_cubic.m.show_in_texture1D(TEXTURE_1D_HSV);
                    computations_cubic.m.show_out_texture1D(TEXTURE_1D_HSV);


                }
            }
        }

        if(gui.but.isChecked())
        {
            computations_cubic.f=get_scalar_field(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),0,1);
            computations_cubic.V=compute_field(computations_cubic.m,computations_cubic.f,gui.method.currentIndex());
            computations_cubic.V_norm=computations_cubic.V;
            computations_cubic.V_norm.normalize();
            gui.canvas.push_obj(&computations_cubic.V_norm,false);


        }
        if(gui.ground_truth.isChecked())
        {
            computations_cubic.GT=compute_ground_truth(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
            computations_cubic.GT_norm=computations_cubic.GT;
            computations_cubic.GT_norm.normalize();
            gui.canvas.push_obj(&computations_cubic.GT_norm,false);

        }


        gui.canvas.updateGL();

    });
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    QSlider::connect(&gui.sl_error_pos, &QSlider::valueChanged, [&]()
    {


        if (gui.Err.isChecked())
        {
            gui.canvas.pop(&computations_cubic.dual_m);
            gui.canvas.pop(&computations_cubic.m_grid);
            if(gui.InsideError.isChecked())
            {
                computations_cubic.f=get_scalar_field(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),0,1);
                computations_cubic.V=compute_field(computations_cubic.m,computations_cubic.f,gui.method.currentIndex());


                bring_the_field_inside(computations_cubic.m,computations_cubic.m_grid,computations_cubic.V,computations_cubic.V_grid,gui.method.currentIndex());
                computations_cubic.GT=compute_ground_truth(computations_cubic.m_grid,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());

                computations_cubic.err=estimate_error(computations_cubic.GT,computations_cubic.V_grid,computations_cubic.m_grid,gui.mode.currentIndex(),LSDD,gui.ErrType.currentIndex(),gui.type_of_vertices.currentIndex());
                computations_cubic.max=find_max_norm(computations_cubic.GT);
                computations_cubic.err_norm=heat_map_normalization(computations_cubic.err,0,computations_cubic.max,gui.sl_error_neg.value(),gui.sl_error_pos.value());
                computations_cubic.err_norm.copy_to_mesh(computations_cubic.m_grid);

                gui.canvas.pop(&computations_cubic.m);
                gui.canvas.push_obj(&computations_cubic.m_grid,false);
                computations_cubic.m_grid.show_in_wireframe(false);
                computations_cubic.m_grid.show_out_wireframe(false);
                computations_cubic.m_grid.show_in_texture1D(TEXTURE_1D_HSV);
                computations_cubic.m_grid.show_out_texture1D(TEXTURE_1D_HSV);
                gui.canvas.push_obj(&computations_cubic.m,false);

            }else{

                computations_cubic.f=get_scalar_field(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),0,1);
                computations_cubic.GT=compute_ground_truth(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
                computations_cubic.V=compute_field(computations_cubic.m,computations_cubic.f,gui.method.currentIndex());
                computations_cubic.err=estimate_error(computations_cubic.GT,computations_cubic.V,computations_cubic.m,gui.mode.currentIndex(),gui.method.currentIndex(),gui.ErrType.currentIndex(),gui.type_of_vertices.currentIndex());
                computations_cubic.max=find_max_norm(computations_cubic.GT);
                computations_cubic.err_norm=heat_map_normalization(computations_cubic.err,0,computations_cubic.max,gui.sl_error_neg.value(),gui.sl_error_pos.value());

                if(gui.method.currentIndex()==0)
                {


                    computations_cubic.err_norm.copy_to_mesh(computations_cubic.dual_m);

                    gui.canvas.pop(&computations_cubic.m);
                    gui.canvas.push_obj(&computations_cubic.dual_m,false);
                    computations_cubic.dual_m.show_in_wireframe(false);
                    computations_cubic.dual_m.show_out_wireframe(false);


                    computations_cubic.dual_m.show_in_texture1D(TEXTURE_1D_HSV);
                    computations_cubic.dual_m.show_out_texture1D(TEXTURE_1D_HSV);
                    gui.canvas.push_obj(&computations_cubic.m,false);
                }else
                {
                    computations_cubic.err_norm.copy_to_mesh(computations_cubic.m);
                    computations_cubic.m.show_in_texture1D(TEXTURE_1D_HSV);
                    computations_cubic.m.show_out_texture1D(TEXTURE_1D_HSV);


                }
            }
        }

        if(gui.but.isChecked())
        {
            computations_cubic.f=get_scalar_field(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),0,1);
            computations_cubic.V=compute_field(computations_cubic.m,computations_cubic.f,gui.method.currentIndex());
            computations_cubic.V_norm=computations_cubic.V;
            computations_cubic.V_norm.normalize();
            gui.canvas.push_obj(&computations_cubic.V_norm,false);


        }
        if(gui.ground_truth.isChecked())
        {
            computations_cubic.GT=compute_ground_truth(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
            computations_cubic.GT_norm=computations_cubic.GT;
            computations_cubic.GT_norm.normalize();
            gui.canvas.push_obj(&computations_cubic.GT_norm,false);

        }


        gui.canvas.updateGL();

    });

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    QCheckBox::connect(&gui.ground_truth, &QCheckBox::stateChanged, [&]()
    {
        if(gui.ground_truth.isChecked())
        {
            computations_cubic.GT=compute_ground_truth(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
            computations_cubic.GT_norm=computations_cubic.GT;
            computations_cubic.GT_norm.normalize();
            gui.canvas.push_obj(&computations_cubic.GT_norm,false);

        }else
        {
            gui.canvas.pop_all_occurrences_of(DRAWABLE_VECTOR_FIELD);

            if(gui.but.isChecked())
            {
                //                computations_cubic.f=get_scalar_field(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),0,0);
                //                computations_cubic.V=compute_field(computations_cubic.m,computations_cubic.f,gui.method.currentIndex());
                //                computations_cubic.V_norm=computations_cubic.V;
                //                computations_cubic.V_norm.normalize();
                gui.canvas.push_obj(&computations_cubic.V_norm,false);


            }
        }

        gui.canvas.updateGL();
    });
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    QCheckBox::connect(&gui.wireframe, &QCheckBox::stateChanged, [&]()
    {
        if(gui.wireframe.isChecked())
        {
            computations_cubic.m.show_in_wireframe(true);
            computations_cubic.m.show_out_wireframe(true);
        }else{
            computations_cubic.m.show_in_wireframe(false);
            computations_cubic.m.show_out_wireframe(false);
        }
        gui.canvas.updateGL();
    });
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    QComboBox::connect(&gui.cb_function, static_cast<void(QComboBox::*)(int)>(&QComboBox::activated), [&]()
    {


        if(gui.but.isChecked())
        {
            computations_cubic.f=get_scalar_field(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),0,0);
            computations_cubic.V=compute_field(computations_cubic.m,computations_cubic.f,gui.method.currentIndex());
            computations_cubic.V_norm=computations_cubic.V;
            computations_cubic.V_norm.normalize();
            gui.canvas.push_obj(&computations_cubic.V_norm,false);


        }
        if(gui.ground_truth.isChecked())
        {
            computations_cubic.GT=compute_ground_truth(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
            computations_cubic.GT_norm=computations_cubic.GT;
            computations_cubic.GT_norm.normalize();
            gui.canvas.push_obj(&computations_cubic.GT_norm,false);

        }

        if (gui.show_heat_map.isChecked())
        {
            if (gui.scalar_field.isChecked())
            {
                gui.canvas.pop(&computations_cubic.dual_m);
                gui.canvas.pop(&computations_cubic.m_grid);



                computations_cubic.f=get_scalar_field(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),0,0);
                find_max_min_values(computations_cubic.f,computations_cubic.max,computations_cubic.min);
                computations_cubic.f_norm=heat_map_normalization(computations_cubic.f,computations_cubic.min,computations_cubic.max);
                computations_cubic.f_norm.copy_to_mesh(computations_cubic.m);
                computations_cubic.m.show_in_texture1D(TEXTURE_1D_HSV);


            }

            if (gui.Err.isChecked())
            {

                gui.canvas.pop(&computations_cubic.dual_m);
                gui.canvas.pop(&computations_cubic.m_grid);
                //computations_cubic.dual_m.clear();

                if(gui.InsideError.isChecked())
                {

                    computations_cubic.f=get_scalar_field(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),0,0);
                    computations_cubic.V=compute_field(computations_cubic.m,computations_cubic.f,gui.method.currentIndex());


                    bring_the_field_inside(computations_cubic.m,computations_cubic.m_grid,computations_cubic.V,computations_cubic.V_grid,gui.method.currentIndex());
                    computations_cubic.GT=compute_ground_truth(computations_cubic.m_grid,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),LSDD);
                    computations_cubic.err=estimate_error(computations_cubic.GT,computations_cubic.V_grid,computations_cubic.m_grid,gui.mode.currentIndex(),LSDD,gui.ErrType.currentIndex(),gui.type_of_vertices.currentIndex());
                    computations_cubic.max=find_max_norm(computations_cubic.GT);
                    computations_cubic.err_norm=heat_map_normalization(computations_cubic.err,0,computations_cubic.max,gui.sl_error_neg.value(),gui.sl_error_pos.value());
                    computations_cubic.err_norm.copy_to_mesh(computations_cubic.m_grid);

                    gui.canvas.pop(&computations_cubic.m);
                    gui.canvas.push_obj(&computations_cubic.m_grid,false);
                    computations_cubic.m_grid.show_in_wireframe(false);
                    computations_cubic.m_grid.show_out_wireframe(false);
                    computations_cubic.m_grid.show_in_texture1D(TEXTURE_1D_HSV);
                    computations_cubic.m.show_in_vert_color();
                    computations_cubic.m.show_out_vert_color();
                    gui.canvas.push_obj(&computations_cubic.m,false);

                }else
                { computations_cubic.f=get_scalar_field(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),0,0);
                    computations_cubic.GT=compute_ground_truth(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
                    computations_cubic.V=compute_field(computations_cubic.m,computations_cubic.f,gui.method.currentIndex());
                    computations_cubic.err=estimate_error(computations_cubic.GT,computations_cubic.V,computations_cubic.m,gui.mode.currentIndex(),gui.method.currentIndex(),gui.ErrType.currentIndex(),gui.type_of_vertices.currentIndex());
                    computations_cubic.max=find_max_norm(computations_cubic.GT);
                    computations_cubic.err_norm=heat_map_normalization(computations_cubic.err,0,computations_cubic.max,gui.sl_error_neg.value(),gui.sl_error_pos.value());

                    if(gui.method.currentIndex()==0)
                    {


                        computations_cubic.err_norm.copy_to_mesh(computations_cubic.dual_m);

                        gui.canvas.pop(&computations_cubic.m);
                        gui.canvas.push_obj(&computations_cubic.dual_m,false);
                        computations_cubic.dual_m.show_in_wireframe(false);
                        computations_cubic.dual_m.show_out_wireframe(false);

                        computations_cubic.dual_m.show_in_texture1D(TEXTURE_1D_HSV);
                        gui.canvas.push_obj(&computations_cubic.m,false);
                    }else
                    {
                        computations_cubic.err_norm.copy_to_mesh(computations_cubic.m);
                        computations_cubic.m.show_in_texture1D(TEXTURE_1D_HSV);


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
            //computations_cubic.f=get_scalar_field(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),0,0);
            computations_cubic.V=compute_field(computations_cubic.m,computations_cubic.f,gui.method.currentIndex());
            computations_cubic.V_norm=computations_cubic.V;
            computations_cubic.V_norm.normalize();
            gui.canvas.push_obj(&computations_cubic.V_norm,false);


        }

        if(gui.ground_truth.isChecked())
        {
            computations_cubic.GT=compute_ground_truth(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
            computations_cubic.GT_norm=computations_cubic.GT;
            computations_cubic.GT_norm.normalize();
            gui.canvas.push_obj(&computations_cubic.GT_norm,false);

        }

        if (gui.show_heat_map.isChecked())
        {

            if (gui.scalar_field.isChecked())
            {
                gui.canvas.pop(&computations_cubic.dual_m);
                gui.canvas.pop(&computations_cubic.m_grid);



                computations_cubic.f=get_scalar_field(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),0,0);
                find_max_min_values(computations_cubic.f,computations_cubic.max,computations_cubic.min);
                computations_cubic.f_norm=heat_map_normalization(computations_cubic.f,computations_cubic.min,computations_cubic.max);
                computations_cubic.f_norm.copy_to_mesh(computations_cubic.m);
                computations_cubic.m.show_in_texture1D(TEXTURE_1D_HSV);


            }
            if (gui.Err.isChecked())
            {

                gui.canvas.pop(&computations_cubic.dual_m);
                gui.canvas.pop(&computations_cubic.m_grid);
                //computations_cubic.dual_m.clear();

                if(gui.InsideError.isChecked())
                {

                    computations_cubic.f=get_scalar_field(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),0,0);
                    computations_cubic.V=compute_field(computations_cubic.m,computations_cubic.f,gui.method.currentIndex());

                    bring_the_field_inside(computations_cubic.m,computations_cubic.m_grid,computations_cubic.V,computations_cubic.V_grid,gui.method.currentIndex());
                    computations_cubic.GT=compute_ground_truth(computations_cubic.m_grid,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),LSDD);
                    computations_cubic.err=estimate_error(computations_cubic.GT,computations_cubic.V_grid,computations_cubic.m_grid,gui.mode.currentIndex(),LSDD,gui.ErrType.currentIndex(),gui.type_of_vertices.currentIndex());
                    computations_cubic.max=find_max_norm(computations_cubic.GT);
                    computations_cubic.err_norm=heat_map_normalization(computations_cubic.err,0,computations_cubic.max,gui.sl_error_neg.value(),gui.sl_error_pos.value());
                    computations_cubic.err_norm.copy_to_mesh(computations_cubic.m_grid);

                    gui.canvas.pop(&computations_cubic.m);
                    gui.canvas.push_obj(&computations_cubic.m_grid,false);
                    computations_cubic.m_grid.show_in_wireframe(false);
                    computations_cubic.m_grid.show_out_wireframe(false);
                    computations_cubic.m_grid.show_in_texture1D(TEXTURE_1D_HSV);
                    computations_cubic.m.show_in_vert_color();
                    computations_cubic.m.show_out_vert_color();
                    gui.canvas.push_obj(&computations_cubic.m,false);


                }else{

                    computations_cubic.f=get_scalar_field(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),0,0);
                    computations_cubic.GT=compute_ground_truth(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
                    computations_cubic.V=compute_field(computations_cubic.m,computations_cubic.f,gui.method.currentIndex());
                    computations_cubic.err=estimate_error(computations_cubic.GT,computations_cubic.V,computations_cubic.m,gui.mode.currentIndex(),gui.method.currentIndex(),gui.ErrType.currentIndex(),gui.type_of_vertices.currentIndex());
                    computations_cubic.max=find_max_norm(computations_cubic.GT);
                    computations_cubic.err_norm=heat_map_normalization(computations_cubic.err,0,computations_cubic.max,gui.sl_error_neg.value(),gui.sl_error_pos.value());

                    if(gui.method.currentIndex()==0)
                    {


                        computations_cubic.err_norm.copy_to_mesh(computations_cubic.dual_m);


                        gui.canvas.pop(&computations_cubic.m);
                        gui.canvas.push_obj(&computations_cubic.dual_m,false);
                        computations_cubic.dual_m.show_in_wireframe(false);
                        computations_cubic.dual_m.show_out_wireframe(false);

                        computations_cubic.dual_m.show_in_texture1D(TEXTURE_1D_HSV);
                        gui.canvas.push_obj(&computations_cubic.m,false);
                    }else
                    {
                        computations_cubic.err_norm.copy_to_mesh(computations_cubic.m);
                        computations_cubic.m.show_in_texture1D(TEXTURE_1D_HSV);


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
        make_tet_mesh(computations_cubic.m,gui.N.value(),gui.choose_tri.currentIndex(),gui.anis.value());
        //dual_mesh(computations_cubic.m, computations_cubic.dual_m,true);
        computations_cubic.dual_m.updateGL();
        computations_cubic.m.slice(gui.s);
        if(gui.method.currentIndex()==0)
        {
            computations_cubic.dual_m.slice(gui.s);
        }


        if(gui.but.isChecked())
        {
            computations_cubic.f=get_scalar_field(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),0,0);
            computations_cubic.V=compute_field(computations_cubic.m,computations_cubic.f,gui.method.currentIndex());
            computations_cubic.V_norm=computations_cubic.V;
            computations_cubic.V_norm.normalize();
            gui.canvas.push_obj(&computations_cubic.V_norm,false);


        }
        if(gui.ground_truth.isChecked())
        {
            computations_cubic.GT=compute_ground_truth(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
            computations_cubic.GT_norm=computations_cubic.GT;
            computations_cubic.GT_norm.normalize();
            gui.canvas.push_obj(&computations_cubic.GT_norm,false);

        }

        if (gui.show_heat_map.isChecked())
        {

            if (gui.scalar_field.isChecked())
            {
                gui.canvas.pop(&computations_cubic.dual_m);
                gui.canvas.pop(&computations_cubic.m_grid);
                // computations_cubic.dual_m.clear();


                computations_cubic.f=get_scalar_field(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),0,0);
                find_max_min_values(computations_cubic.f,computations_cubic.max,computations_cubic.min);
                computations_cubic.f_norm=heat_map_normalization(computations_cubic.f,computations_cubic.min,computations_cubic.max);
                computations_cubic.f_norm.copy_to_mesh(computations_cubic.m);
                computations_cubic.m.show_in_texture1D(TEXTURE_1D_HSV);


            }
            if (gui.Err.isChecked())
            {
                gui.canvas.pop(&computations_cubic.dual_m);
                gui.canvas.pop(&computations_cubic.m_grid);
                //computations_cubic.dual_m.clear();
                if(gui.InsideError.isChecked())

                {
                    computations_cubic.f=get_scalar_field(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),0,0);
                    computations_cubic.V=compute_field(computations_cubic.m,computations_cubic.f,gui.method.currentIndex());

                    bring_the_field_inside(computations_cubic.m,computations_cubic.m_grid,computations_cubic.V,computations_cubic.V_grid,gui.method.currentIndex());
                    computations_cubic.GT=compute_ground_truth(computations_cubic.m_grid,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),LSDD);
                    computations_cubic.err=estimate_error(computations_cubic.GT,computations_cubic.V_grid,computations_cubic.m_grid,gui.mode.currentIndex(),LSDD,gui.ErrType.currentIndex(),gui.type_of_vertices.currentIndex());
                    computations_cubic.max=find_max_norm(computations_cubic.GT);
                    computations_cubic.err_norm=heat_map_normalization(computations_cubic.err,0,computations_cubic.max,gui.sl_error_neg.value(),gui.sl_error_pos.value());
                    computations_cubic.err_norm.copy_to_mesh(computations_cubic.m_grid);

                    gui.canvas.pop(&computations_cubic.m);
                    gui.canvas.push_obj(&computations_cubic.m_grid,false);
                    computations_cubic.m_grid.show_in_wireframe(false);
                    computations_cubic.m_grid.show_out_wireframe(false);
                    computations_cubic.m_grid.show_in_texture1D(TEXTURE_1D_HSV);
                    computations_cubic.m.show_in_vert_color();
                    computations_cubic.m.show_out_vert_color();
                    gui.canvas.push_obj(&computations_cubic.m,false);
                }
                else{
                    computations_cubic.f=get_scalar_field(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),0,0);
                    computations_cubic.GT=compute_ground_truth(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
                    computations_cubic.V=compute_field(computations_cubic.m,computations_cubic.f,gui.method.currentIndex());
                    computations_cubic.err=estimate_error(computations_cubic.GT,computations_cubic.V,computations_cubic.m,gui.mode.currentIndex(),gui.method.currentIndex(),gui.ErrType.currentIndex(),gui.type_of_vertices.currentIndex());
                    computations_cubic.max=find_max_norm(computations_cubic.GT);
                    computations_cubic.err_norm=heat_map_normalization(computations_cubic.err,0,computations_cubic.max,gui.sl_error_neg.value(),gui.sl_error_pos.value());

                    if(gui.method.currentIndex()==0)
                    {


                        computations_cubic.err_norm.copy_to_mesh(computations_cubic.dual_m);

                        gui.canvas.pop(&computations_cubic.m);
                        gui.canvas.push_obj(&computations_cubic.dual_m,false);
                        computations_cubic.dual_m.show_in_wireframe(false);
                        computations_cubic.dual_m.show_out_wireframe(false);

                        computations_cubic.dual_m.show_in_texture1D(TEXTURE_1D_HSV);
                        gui.canvas.push_obj(&computations_cubic.m,false);
                    }else
                    {
                        computations_cubic.err_norm.copy_to_mesh(computations_cubic.m);
                        computations_cubic.m.show_in_texture1D(TEXTURE_1D_HSV);


                    }
                }
            }
        }




        gui.canvas.updateGL();

    });
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    QSpinBox::connect(&gui.N,static_cast<void(QSpinBox::*)(int)>(&QSpinBox::valueChanged), [&]()
    {


        make_tet_mesh(computations_cubic.m,gui.N.value(),gui.choose_tri.currentIndex());
        //dual_mesh(computations_cubic.m, computations_cubic.dual_m,true);
        computations_cubic.dual_m.updateGL();
        computations_cubic.m.slice(gui.s);
        if(gui.method.currentIndex()==0)
        {
            computations_cubic.dual_m.slice(gui.s);
        }




        if(gui.but.isChecked())
        {
            computations_cubic.f=get_scalar_field(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),0,0);
            computations_cubic.V=compute_field(computations_cubic.m,computations_cubic.f,gui.method.currentIndex());
            computations_cubic.V_norm=computations_cubic.V;
            computations_cubic.V_norm.normalize();
            gui.canvas.push_obj(&computations_cubic.V_norm,false);


        }

        if(gui.ground_truth.isChecked())
        {
            computations_cubic.GT=compute_ground_truth(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
            computations_cubic.GT_norm=computations_cubic.GT;
            computations_cubic.GT_norm.normalize();
            gui.canvas.push_obj(&computations_cubic.GT_norm,false);

        }
        if (gui.show_heat_map.isChecked())
        {

            if (gui.scalar_field.isChecked())
            {


                computations_cubic.f=get_scalar_field(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),0,0);
                find_max_min_values(computations_cubic.f,computations_cubic.max,computations_cubic.min);
                computations_cubic.f_norm=heat_map_normalization(computations_cubic.f,computations_cubic.min,computations_cubic.max);
                computations_cubic.f_norm.copy_to_mesh(computations_cubic.m);
                computations_cubic.m.show_in_texture1D(TEXTURE_1D_HSV);


            }
            if (gui.Err.isChecked())
            {

                gui.canvas.pop(&computations_cubic.dual_m);
                gui.canvas.pop(&computations_cubic.m_grid);
                //computations_cubic.dual_m.clear();

                if(gui.InsideError.isChecked())
                {

                    computations_cubic.f=get_scalar_field(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),0,0);
                    computations_cubic.V=compute_field(computations_cubic.m,computations_cubic.f,gui.method.currentIndex());

                    bring_the_field_inside(computations_cubic.m,computations_cubic.m_grid,computations_cubic.V,computations_cubic.V_grid,gui.method.currentIndex());
                    computations_cubic.GT=compute_ground_truth(computations_cubic.m_grid,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),LSDD);
                    computations_cubic.err=estimate_error(computations_cubic.GT,computations_cubic.V_grid,computations_cubic.m_grid,gui.mode.currentIndex(),LSDD,gui.ErrType.currentIndex(),gui.type_of_vertices.currentIndex());
                    computations_cubic.max=find_max_norm(computations_cubic.GT);
                    computations_cubic.err_norm=heat_map_normalization(computations_cubic.err,0,computations_cubic.max,gui.sl_error_neg.value(),gui.sl_error_pos.value());
                    computations_cubic.err_norm.copy_to_mesh(computations_cubic.m_grid);

                    gui.canvas.pop(&computations_cubic.m);
                    gui.canvas.push_obj(&computations_cubic.m_grid,false);
                    computations_cubic.m_grid.show_in_wireframe(false);
                    computations_cubic.m_grid.show_out_wireframe(false);
                    computations_cubic.m_grid.show_in_texture1D(TEXTURE_1D_HSV);
                    computations_cubic.m.show_in_vert_color();
                    computations_cubic.m.show_out_vert_color();
                    gui.canvas.push_obj(&computations_cubic.m,false);

                }else{
                    computations_cubic.f=get_scalar_field(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),0,0);
                    computations_cubic.GT=compute_ground_truth(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
                    computations_cubic.V=compute_field(computations_cubic.m,computations_cubic.f,gui.method.currentIndex());
                    computations_cubic.err=estimate_error(computations_cubic.GT,computations_cubic.V,computations_cubic.m,gui.mode.currentIndex(),gui.method.currentIndex(),gui.ErrType.currentIndex(),gui.type_of_vertices.currentIndex());
                    computations_cubic.max=find_max_norm(computations_cubic.GT);
                    computations_cubic.err_norm=heat_map_normalization(computations_cubic.err,0,computations_cubic.max,gui.sl_error_neg.value(),gui.sl_error_pos.value());

                    if(gui.method.currentIndex()==0)
                    {


                        computations_cubic.err_norm.copy_to_mesh(computations_cubic.dual_m);

                        gui.canvas.pop(&computations_cubic.m);
                        gui.canvas.push_obj(&computations_cubic.dual_m,false);
                        computations_cubic.dual_m.show_in_wireframe(false);
                        computations_cubic.dual_m.show_out_wireframe(false);

                        computations_cubic.dual_m.show_in_texture1D(TEXTURE_1D_HSV);
                        gui.canvas.push_obj(&computations_cubic.m,false);
                    }else
                    {
                        computations_cubic.err_norm.copy_to_mesh(computations_cubic.m);
                        computations_cubic.m.show_in_texture1D(TEXTURE_1D_HSV);


                    }
                }
            }
        }
        gui.canvas.updateGL();

    });
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    QSpinBox::connect(&gui.anis,static_cast<void(QSpinBox::*)(int)>(&QSpinBox::valueChanged), [&]()
    {
        make_tet_mesh(computations_cubic.m,gui.N.value(),gui.choose_tri.currentIndex(),gui.anis.value());


        if(gui.but.isChecked())
        {
            computations_cubic.f=get_scalar_field(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),0,0);
            computations_cubic.V=compute_field(computations_cubic.m,computations_cubic.f,gui.method.currentIndex());
            computations_cubic.V_norm=computations_cubic.V;
            computations_cubic.V_norm.normalize();
            gui.canvas.push_obj(&computations_cubic.V_norm,false);


        }

        if(gui.ground_truth.isChecked())
        {
            computations_cubic.GT=compute_ground_truth(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
            computations_cubic.GT_norm=computations_cubic.GT;
            computations_cubic.GT_norm.normalize();
            gui.canvas.push_obj(&computations_cubic.GT_norm,false);

        }
        if (gui.show_heat_map.isChecked())
        {

            if (gui.scalar_field.isChecked())
            {
                gui.canvas.pop(&computations_cubic.dual_m);
                gui.canvas.pop(&computations_cubic.m_grid);
                //computations_cubic.dual_m.clear();


                computations_cubic.f=get_scalar_field(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),0,0);
                find_max_min_values(computations_cubic.f,computations_cubic.max,computations_cubic.min);
                computations_cubic.f_norm=heat_map_normalization(computations_cubic.f,computations_cubic.min,computations_cubic.max);
                computations_cubic.f_norm.copy_to_mesh(computations_cubic.m);
                computations_cubic.m.show_in_texture1D(TEXTURE_1D_HSV);


            }
            if (gui.Err.isChecked())
            {

                gui.canvas.pop(&computations_cubic.dual_m);
                gui.canvas.pop(&computations_cubic.m_grid);
                //computations_cubic.dual_m.clear();

                if(gui.InsideError.isChecked())
                {

                    computations_cubic.f=get_scalar_field(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),0,0);
                    computations_cubic.V=compute_field(computations_cubic.m,computations_cubic.f,gui.method.currentIndex());

                    bring_the_field_inside(computations_cubic.m,computations_cubic.m_grid,computations_cubic.V,computations_cubic.V_grid,gui.method.currentIndex());
                    computations_cubic.GT=compute_ground_truth(computations_cubic.m_grid,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),LSDD);
                    computations_cubic.err=estimate_error(computations_cubic.GT,computations_cubic.V_grid,computations_cubic.m_grid,gui.mode.currentIndex(),LSDD,gui.ErrType.currentIndex(),gui.type_of_vertices.currentIndex());
                    computations_cubic.max=find_max_norm(computations_cubic.GT);
                    computations_cubic.err_norm=heat_map_normalization(computations_cubic.err,0,computations_cubic.max,gui.sl_error_neg.value(),gui.sl_error_pos.value());
                    computations_cubic.err_norm.copy_to_mesh(computations_cubic.m_grid);

                    gui.canvas.pop(&computations_cubic.m);
                    gui.canvas.push_obj(&computations_cubic.m_grid,false);
                    computations_cubic.m_grid.show_in_wireframe(false);
                    computations_cubic.m_grid.show_out_wireframe(false);
                    computations_cubic.m_grid.show_in_texture1D(TEXTURE_1D_HSV);
                    computations_cubic.m.show_in_vert_color();
                    computations_cubic.m.show_out_vert_color();
                    gui.canvas.push_obj(&computations_cubic.m,false);

                }else{
                    computations_cubic.f=get_scalar_field(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),0,0);
                    computations_cubic.GT=compute_ground_truth(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
                    computations_cubic.V=compute_field(computations_cubic.m,computations_cubic.f,gui.method.currentIndex());
                    computations_cubic.err=estimate_error(computations_cubic.GT,computations_cubic.V,computations_cubic.m,gui.mode.currentIndex(),gui.method.currentIndex(),gui.ErrType.currentIndex(),gui.type_of_vertices.currentIndex());
                    computations_cubic.max=find_max_norm(computations_cubic.GT);
                    computations_cubic.err_norm=heat_map_normalization(computations_cubic.err,0,computations_cubic.max,gui.sl_error_neg.value(),gui.sl_error_pos.value());

                    if(gui.method.currentIndex()==0)
                    {


                        computations_cubic.err_norm.copy_to_mesh(computations_cubic.dual_m);

                        gui.canvas.pop(&computations_cubic.m);
                        gui.canvas.push_obj(&computations_cubic.dual_m,false);
                        computations_cubic.dual_m.show_in_wireframe(false);
                        computations_cubic.dual_m.show_out_wireframe(false);

                        computations_cubic.dual_m.show_in_texture1D(TEXTURE_1D_HSV);
                        gui.canvas.push_obj(&computations_cubic.m,false);
                    }else
                    {
                        computations_cubic.err_norm.copy_to_mesh(computations_cubic.m);
                        computations_cubic.m.show_in_texture1D(TEXTURE_1D_HSV);


                    }
                }
            }
        }
        gui.canvas.updateGL();
    });


    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    QSpinBox::connect(&gui.a,static_cast<void(QSpinBox::*)(int)>(&QSpinBox::valueChanged), [&]()
    {
        computations_cubic.f=get_scalar_field(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),0,0);
        computations_cubic.GT=compute_ground_truth(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());

        if(gui.but.isChecked())
        {
            computations_cubic.f=get_scalar_field(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),0,0);
            computations_cubic.V=compute_field(computations_cubic.m,computations_cubic.f,gui.method.currentIndex());
            computations_cubic.V_norm=computations_cubic.V;
            computations_cubic.V_norm.normalize();
            gui.canvas.push_obj(&computations_cubic.V_norm,false);


        }
        if(gui.ground_truth.isChecked())
        {
            computations_cubic.GT=compute_ground_truth(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
            computations_cubic.GT_norm=computations_cubic.GT;
            computations_cubic.GT_norm.normalize();
            gui.canvas.push_obj(&computations_cubic.GT_norm,false);

        }

        if (gui.show_heat_map.isChecked())
        {
            if (gui.scalar_field.isChecked())
            {
                gui.canvas.pop(&computations_cubic.dual_m);
                gui.canvas.pop(&computations_cubic.m_grid);
                //computations_cubic.dual_m.clear();


                computations_cubic.f=get_scalar_field(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),0,0);
                find_max_min_values(computations_cubic.f,computations_cubic.max,computations_cubic.min);
                computations_cubic.f_norm=heat_map_normalization(computations_cubic.f,computations_cubic.min,computations_cubic.max);
                computations_cubic.f_norm.copy_to_mesh(computations_cubic.m);
                computations_cubic.m.show_in_texture1D(TEXTURE_1D_HSV);


            }
            if (gui.Err.isChecked())
            {

                gui.canvas.pop(&computations_cubic.dual_m);
                gui.canvas.pop(&computations_cubic.m_grid);
                //computations_cubic.dual_m.clear();

                if(gui.InsideError.isChecked())
                {

                    computations_cubic.f=get_scalar_field(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),0,0);
                    computations_cubic.V=compute_field(computations_cubic.m,computations_cubic.f,gui.method.currentIndex());

                    bring_the_field_inside(computations_cubic.m,computations_cubic.m_grid,computations_cubic.V,computations_cubic.V_grid,gui.method.currentIndex());
                    computations_cubic.GT=compute_ground_truth(computations_cubic.m_grid,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),LSDD);
                    computations_cubic.err=estimate_error(computations_cubic.GT,computations_cubic.V_grid,computations_cubic.m_grid,gui.mode.currentIndex(),LSDD,gui.ErrType.currentIndex(),gui.type_of_vertices.currentIndex());
                    computations_cubic.max=find_max_norm(computations_cubic.GT);
                    computations_cubic.err_norm=heat_map_normalization(computations_cubic.err,0,computations_cubic.max,gui.sl_error_neg.value(),gui.sl_error_pos.value());
                    computations_cubic.err_norm.copy_to_mesh(computations_cubic.m_grid);

                    gui.canvas.pop(&computations_cubic.m);
                    gui.canvas.push_obj(&computations_cubic.m_grid,false);
                    computations_cubic.m_grid.show_in_wireframe(false);
                    computations_cubic.m_grid.show_out_wireframe(false);
                    computations_cubic.m_grid.show_in_texture1D(TEXTURE_1D_HSV);
                    computations_cubic.m.show_in_vert_color();
                    computations_cubic.m.show_out_vert_color();
                    gui.canvas.push_obj(&computations_cubic.m,false);
                }else{
                    computations_cubic.f=get_scalar_field(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),0,0);
                    computations_cubic.GT=compute_ground_truth(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
                    computations_cubic.V=compute_field(computations_cubic.m,computations_cubic.f,gui.method.currentIndex());
                    computations_cubic.err=estimate_error(computations_cubic.GT,computations_cubic.V,computations_cubic.m,gui.mode.currentIndex(),gui.method.currentIndex(),gui.ErrType.currentIndex(),gui.type_of_vertices.currentIndex());
                    computations_cubic.max=find_max_norm(computations_cubic.GT);
                    computations_cubic.err_norm=heat_map_normalization(computations_cubic.err,0,computations_cubic.max,gui.sl_error_neg.value(),gui.sl_error_pos.value());

                    if(gui.method.currentIndex()==0)
                    {


                        computations_cubic.err_norm.copy_to_mesh(computations_cubic.dual_m);

                        gui.canvas.pop(&computations_cubic.m);
                        gui.canvas.push_obj(&computations_cubic.dual_m,false);
                        computations_cubic.dual_m.show_in_wireframe(false);
                        computations_cubic.dual_m.show_out_wireframe(false);

                        computations_cubic.dual_m.show_in_texture1D(TEXTURE_1D_HSV);
                        gui.canvas.push_obj(&computations_cubic.m,false);
                    }else
                    {
                        computations_cubic.err_norm.copy_to_mesh(computations_cubic.m);
                        computations_cubic.m.show_in_texture1D(TEXTURE_1D_HSV);


                    }
                }}

        }
        gui.canvas.updateGL();




    });
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    QSpinBox::connect(&gui.b,static_cast<void(QSpinBox::*)(int)>(&QSpinBox::valueChanged), [&]()
    {
        computations_cubic.f=get_scalar_field(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),0,0);
        computations_cubic.GT=compute_ground_truth(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());

        if(gui.but.isChecked())
        {
            computations_cubic.f=get_scalar_field(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),0,0);
            computations_cubic.V=compute_field(computations_cubic.m,computations_cubic.f,gui.method.currentIndex());
            computations_cubic.V_norm=computations_cubic.V;
            computations_cubic.V_norm.normalize();
            gui.canvas.push_obj(&computations_cubic.V_norm,false);


        }
        if(gui.ground_truth.isChecked())
        {
            computations_cubic.GT=compute_ground_truth(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
            computations_cubic.GT_norm=computations_cubic.GT;
            computations_cubic.GT_norm.normalize();
            gui.canvas.push_obj(&computations_cubic.GT_norm,false);

        }
        if (gui.show_heat_map.isChecked())
        {

            if (gui.scalar_field.isChecked())
            {
                gui.canvas.pop(&computations_cubic.dual_m);
                gui.canvas.pop(&computations_cubic.m_grid);
                //computations_cubic.dual_m.clear();


                computations_cubic.f=get_scalar_field(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),0,0);
                find_max_min_values(computations_cubic.f,computations_cubic.max,computations_cubic.min);
                computations_cubic.f_norm=heat_map_normalization(computations_cubic.f,computations_cubic.min,computations_cubic.max);
                computations_cubic.f_norm.copy_to_mesh(computations_cubic.m);
                computations_cubic.m.show_in_texture1D(TEXTURE_1D_HSV);


            }
            if (gui.Err.isChecked())
            {

                gui.canvas.pop(&computations_cubic.dual_m);
                gui.canvas.pop(&computations_cubic.m_grid);
                //computations_cubic.dual_m.clear();

                if(gui.InsideError.isChecked())
                {

                    computations_cubic.f=get_scalar_field(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),0,0);
                    computations_cubic.V=compute_field(computations_cubic.m,computations_cubic.f,gui.method.currentIndex());

                    bring_the_field_inside(computations_cubic.m,computations_cubic.m_grid,computations_cubic.V,computations_cubic.V_grid,gui.method.currentIndex());
                    computations_cubic.GT=compute_ground_truth(computations_cubic.m_grid,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),LSDD);
                    computations_cubic.err=estimate_error(computations_cubic.GT,computations_cubic.V_grid,computations_cubic.m_grid,gui.mode.currentIndex(),LSDD,gui.ErrType.currentIndex(),gui.type_of_vertices.currentIndex());
                    computations_cubic.max=find_max_norm(computations_cubic.GT);
                    computations_cubic.err_norm=heat_map_normalization(computations_cubic.err,0,computations_cubic.max,gui.sl_error_neg.value(),gui.sl_error_pos.value());
                    computations_cubic.err_norm.copy_to_mesh(computations_cubic.m_grid);

                    gui.canvas.pop(&computations_cubic.m);
                    gui.canvas.push_obj(&computations_cubic.m_grid,false);
                    computations_cubic.m_grid.show_in_wireframe(false);
                    computations_cubic.m_grid.show_out_wireframe(false);
                    computations_cubic.m_grid.show_in_texture1D(TEXTURE_1D_HSV);
                    computations_cubic.m.show_in_vert_color();
                    computations_cubic.m.show_out_vert_color();
                    gui.canvas.push_obj(&computations_cubic.m,false);


                }else{

                    computations_cubic.f=get_scalar_field(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),0,0);
                    computations_cubic.GT=compute_ground_truth(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
                    computations_cubic.V=compute_field(computations_cubic.m,computations_cubic.f,gui.method.currentIndex());
                    computations_cubic.err=estimate_error(computations_cubic.GT,computations_cubic.V,computations_cubic.m,gui.mode.currentIndex(),gui.method.currentIndex(),gui.ErrType.currentIndex(),gui.type_of_vertices.currentIndex());
                    computations_cubic.max=find_max_norm(computations_cubic.GT);
                    computations_cubic.err_norm=heat_map_normalization(computations_cubic.err,0,computations_cubic.max,gui.sl_error_neg.value(),gui.sl_error_pos.value());

                    if(gui.method.currentIndex()==0)
                    {


                        computations_cubic.err_norm.copy_to_mesh(computations_cubic.dual_m);

                        gui.canvas.pop(&computations_cubic.m);
                        gui.canvas.push_obj(&computations_cubic.dual_m,false);
                        computations_cubic.dual_m.show_in_wireframe(false);
                        computations_cubic.dual_m.show_out_wireframe(false);

                        computations_cubic.dual_m.show_in_texture1D(TEXTURE_1D_HSV);
                        gui.canvas.push_obj(&computations_cubic.m,false);
                    }else
                    {
                        computations_cubic.err_norm.copy_to_mesh(computations_cubic.m);
                        computations_cubic.m.show_in_texture1D(TEXTURE_1D_HSV);


                    }
                }}
        }

        gui.canvas.updateGL();
    });
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    QSlider::connect(&gui.sl_slice_x, &QSlider::valueChanged, [&]()
    {

        gui.s.X_thresh = float(gui.sl_slice_x.value()) / float(gui.sl_slice_x.maximum());
        gui.s.Y_thresh = float(gui.sl_slice_y.value()) / float(gui.sl_slice_y.maximum());
        gui.s.Z_thresh = float(gui.sl_slice_z.value()) / float(gui.sl_slice_z.maximum());


        gui.s.X_sign   = gui.cb_slice_flip_x.isChecked() ? GEQ : LEQ;
        gui.s.Y_sign   = gui.cb_slice_flip_y.isChecked() ? GEQ : LEQ;
        gui.s.Z_sign   = gui.cb_slice_flip_z.isChecked() ? GEQ : LEQ;

        computations_cubic.m.slice(gui.s);
        if(gui.method.currentIndex()==0)
        {
            computations_cubic.dual_m.slice(gui.s);
        }
        if(gui.InsideError.isChecked())
        {
            computations_cubic.m_grid.slice(gui.s);
        }
        gui.canvas.updateGL();
    });

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    QSlider::connect(&gui.sl_slice_y, &QSlider::valueChanged, [&]()
    {

        gui.s.X_thresh = float(gui.sl_slice_x.value()) / float(gui.sl_slice_x.maximum());
        gui.s.Y_thresh = float(gui.sl_slice_y.value()) / float(gui.sl_slice_y.maximum());
        gui.s.Z_thresh = float(gui.sl_slice_z.value()) / float(gui.sl_slice_z.maximum());


        gui.s.X_sign   = gui.cb_slice_flip_x.isChecked() ? GEQ : LEQ;
        gui.s.Y_sign   = gui.cb_slice_flip_y.isChecked() ? GEQ : LEQ;
        gui.s.Z_sign   = gui.cb_slice_flip_z.isChecked() ? GEQ : LEQ;

        computations_cubic.m.slice(gui.s);
        if(gui.method.currentIndex()==0)
        {
            computations_cubic.dual_m.slice(gui.s);
        }
        if(gui.InsideError.isChecked())
        {
            computations_cubic.m_grid.slice(gui.s);
        }
        gui.canvas.updateGL();
    });

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    QSlider::connect(&gui.sl_slice_z, &QSlider::valueChanged, [&]()
    {

        gui.s.X_thresh = float(gui.sl_slice_x.value()) / float(gui.sl_slice_x.maximum());
        gui.s.Y_thresh = float(gui.sl_slice_y.value()) / float(gui.sl_slice_y.maximum());
        gui.s.Z_thresh = float(gui.sl_slice_z.value()) / float(gui.sl_slice_z.maximum());


        gui.s.X_sign   = gui.cb_slice_flip_x.isChecked() ? GEQ : LEQ;
        gui.s.Y_sign   = gui.cb_slice_flip_y.isChecked() ? GEQ : LEQ;
        gui.s.Z_sign   = gui.cb_slice_flip_z.isChecked() ? GEQ : LEQ;

        computations_cubic.m.slice(gui.s);
        if(gui.method.currentIndex()==0)
        {
            computations_cubic.dual_m.slice(gui.s);
        }
        if(gui.InsideError.isChecked())
        {
            computations_cubic.m_grid.slice(gui.s);
        }
        gui.canvas.updateGL();
    });

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    QCheckBox::connect(&gui.cb_slice_flip_x, &QCheckBox::toggled, [&]()
    {
        gui.s.X_thresh = float(gui.sl_slice_x.value()) / float(gui.sl_slice_x.maximum());
        gui.s.Y_thresh = float(gui.sl_slice_y.value()) / float(gui.sl_slice_y.maximum());
        gui.s.Z_thresh = float(gui.sl_slice_z.value()) / float(gui.sl_slice_z.maximum());


        gui.s.X_sign   = gui.cb_slice_flip_x.isChecked() ? GEQ : LEQ;
        gui.s.Y_sign   = gui.cb_slice_flip_y.isChecked() ? GEQ : LEQ;
        gui.s.Z_sign   = gui.cb_slice_flip_z.isChecked() ? GEQ : LEQ;

        computations_cubic.m.slice(gui.s);
        if(gui.method.currentIndex()==0)
        {
            computations_cubic.dual_m.slice(gui.s);
        }
        if(gui.InsideError.isChecked())
        {
            computations_cubic.m_grid.slice(gui.s);
        }
        gui.canvas.updateGL();
    });

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    QCheckBox::connect(&gui.cb_slice_flip_y, &QCheckBox::toggled, [&]()
    {
        gui.s.X_thresh = float(gui.sl_slice_x.value()) / float(gui.sl_slice_x.maximum());
        gui.s.Y_thresh = float(gui.sl_slice_y.value()) / float(gui.sl_slice_y.maximum());
        gui.s.Z_thresh = float(gui.sl_slice_z.value()) / float(gui.sl_slice_z.maximum());


        gui.s.X_sign   = gui.cb_slice_flip_x.isChecked() ? GEQ : LEQ;
        gui.s.Y_sign   = gui.cb_slice_flip_y.isChecked() ? GEQ : LEQ;
        gui.s.Z_sign   = gui.cb_slice_flip_z.isChecked() ? GEQ : LEQ;

        computations_cubic.m.slice(gui.s);
        if(gui.method.currentIndex()==0)
        {
            computations_cubic.dual_m.slice(gui.s);
        }
        if(gui.InsideError.isChecked())
        {
            computations_cubic.m_grid.slice(gui.s);
        }
        gui.canvas.updateGL();
    });

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    QCheckBox::connect(&gui.cb_slice_flip_z, &QCheckBox::toggled, [&]()
    {
        gui.s.X_thresh = float(gui.sl_slice_x.value()) / float(gui.sl_slice_x.maximum());
        gui.s.Y_thresh = float(gui.sl_slice_y.value()) / float(gui.sl_slice_y.maximum());
        gui.s.Z_thresh = float(gui.sl_slice_z.value()) / float(gui.sl_slice_z.maximum());


        gui.s.X_sign   = gui.cb_slice_flip_x.isChecked() ? GEQ : LEQ;
        gui.s.Y_sign   = gui.cb_slice_flip_y.isChecked() ? GEQ : LEQ;
        gui.s.Z_sign   = gui.cb_slice_flip_z.isChecked() ? GEQ : LEQ;

        computations_cubic.m.slice(gui.s);
        if(gui.method.currentIndex()==0)
        {
            computations_cubic.dual_m.slice(gui.s);
        }
        if(gui.InsideError.isChecked())
        {
            computations_cubic.m_grid.slice(gui.s);
        }
        gui.canvas.updateGL();
    });

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    QSpinBox::connect(&gui.c,static_cast<void(QSpinBox::*)(int)>(&QSpinBox::valueChanged), [&]()
    {
        computations_cubic.f=get_scalar_field(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),0,0);
        computations_cubic.GT=compute_ground_truth(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());

        if(gui.but.isChecked())
        {
            computations_cubic.f=get_scalar_field(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),0,0);
            computations_cubic.V=compute_field(computations_cubic.m,computations_cubic.f,gui.method.currentIndex());
            computations_cubic.V_norm=computations_cubic.V;
            computations_cubic.V_norm.normalize();
            gui.canvas.push_obj(&computations_cubic.V_norm,false);


        }
        if(gui.ground_truth.isChecked())
        {
            computations_cubic.GT=compute_ground_truth(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
            computations_cubic.GT_norm=computations_cubic.GT;
            computations_cubic.GT_norm.normalize();
            gui.canvas.push_obj(&computations_cubic.GT_norm,false);
            gui.canvas.updateGL();
        }
        if (gui.show_heat_map.isChecked())
        {
            if (gui.scalar_field.isChecked())
            {
                gui.canvas.pop(&computations_cubic.dual_m);
                gui.canvas.pop(&computations_cubic.m_grid);
                //computations_cubic.dual_m.clear();


                computations_cubic.f=get_scalar_field(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),0,0);
                find_max_min_values(computations_cubic.f,computations_cubic.max,computations_cubic.min);
                computations_cubic.f_norm=heat_map_normalization(computations_cubic.f,computations_cubic.min,computations_cubic.max);
                computations_cubic.f_norm.copy_to_mesh(computations_cubic.m);
                computations_cubic.m.show_in_texture1D(TEXTURE_1D_HSV);


            }
            if (gui.Err.isChecked())
            {

                gui.canvas.pop(&computations_cubic.dual_m);
                gui.canvas.pop(&computations_cubic.m_grid);
                //computations_cubic.dual_m.clear();

                if(gui.InsideError.isChecked())
                {

                    computations_cubic.f=get_scalar_field(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),0,0);
                    computations_cubic.V=compute_field(computations_cubic.m,computations_cubic.f,gui.method.currentIndex());
                    //computations_cubic.GT_grid=compute_values_on_grid(computations_cubic.m_grid, gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());

                    bring_the_field_inside(computations_cubic.m,computations_cubic.m_grid,computations_cubic.V,computations_cubic.V_grid,gui.method.currentIndex());
                    computations_cubic.GT=compute_ground_truth(computations_cubic.m_grid,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),LSDD);
                    computations_cubic.err=estimate_error(computations_cubic.GT,computations_cubic.V_grid,computations_cubic.m_grid,gui.mode.currentIndex(),LSDD,gui.ErrType.currentIndex(),gui.type_of_vertices.currentIndex());
                    computations_cubic.max=find_max_norm(computations_cubic.GT);
                    computations_cubic.err_norm=heat_map_normalization(computations_cubic.err,0,computations_cubic.max,gui.sl_error_neg.value(),gui.sl_error_pos.value());
                    computations_cubic.err_norm.copy_to_mesh(computations_cubic.m_grid);

                    gui.canvas.pop(&computations_cubic.m);
                    gui.canvas.push_obj(&computations_cubic.m_grid,false);
                    computations_cubic.m_grid.show_in_wireframe(false);
                    computations_cubic.m_grid.show_out_wireframe(false);
                    computations_cubic.m_grid.show_in_texture1D(TEXTURE_1D_HSV);
                    computations_cubic.m.show_in_vert_color();
                    computations_cubic.m.show_out_vert_color();
                    gui.canvas.push_obj(&computations_cubic.m,false);

                }else{
                    computations_cubic.f=get_scalar_field(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),0,0);
                    computations_cubic.GT=compute_ground_truth(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
                    computations_cubic.V=compute_field(computations_cubic.m,computations_cubic.f,gui.method.currentIndex());
                    computations_cubic.err=estimate_error(computations_cubic.GT,computations_cubic.V,computations_cubic.m,gui.mode.currentIndex(),gui.method.currentIndex(),gui.ErrType.currentIndex(),gui.type_of_vertices.currentIndex());
                    computations_cubic.max=find_max_norm(computations_cubic.GT);
                    computations_cubic.err_norm=heat_map_normalization(computations_cubic.err,0,computations_cubic.max,gui.sl_error_neg.value(),gui.sl_error_pos.value());

                    if(gui.method.currentIndex()==0)
                    {


                        computations_cubic.err_norm.copy_to_mesh(computations_cubic.dual_m);

                        gui.canvas.pop(&computations_cubic.m);
                        gui.canvas.push_obj(&computations_cubic.dual_m,false);
                        computations_cubic.dual_m.show_in_wireframe(false);
                        computations_cubic.dual_m.show_out_wireframe(false);

                        computations_cubic.dual_m.show_in_texture1D(TEXTURE_1D_HSV);
                        gui.canvas.push_obj(&computations_cubic.m,false);
                    }else
                    {
                        computations_cubic.err_norm.copy_to_mesh(computations_cubic.m);
                        computations_cubic.m.show_in_texture1D(TEXTURE_1D_HSV);


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
            gui.canvas.pop(&computations_cubic.dual_m);
            gui.canvas.pop(&computations_cubic.m_grid);

            gui.error.setDisabled(true);
            computations_cubic.f=get_scalar_field(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),0,0);
            find_max_min_values(computations_cubic.f,computations_cubic.max,computations_cubic.min);
            computations_cubic.f_norm=heat_map_normalization(computations_cubic.f,computations_cubic.min,computations_cubic.max);

            computations_cubic.f_norm.copy_to_mesh(computations_cubic.m);
            computations_cubic.m.show_in_texture1D(TEXTURE_1D_HSV);
            if(gui.but.isChecked())
            {
                computations_cubic.f=get_scalar_field(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),0,0);
                computations_cubic.V=compute_field(computations_cubic.m,computations_cubic.f,gui.method.currentIndex());
                computations_cubic.V_norm=computations_cubic.V;
                computations_cubic.V_norm.normalize();
                gui.canvas.push_obj(&computations_cubic.V_norm,false);


            }
            if(gui.ground_truth.isChecked())
            {
                computations_cubic.GT=compute_ground_truth(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
                computations_cubic.GT_norm=computations_cubic.GT;
                computations_cubic.GT_norm.normalize();
                gui.canvas.push_obj(&computations_cubic.GT_norm,false);

            }
            if(gui.but.isChecked())
            {
                computations_cubic.f=get_scalar_field(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),0,0);
                computations_cubic.V=compute_field(computations_cubic.m,computations_cubic.f,gui.method.currentIndex());
                computations_cubic.V_norm=computations_cubic.V;
                computations_cubic.V_norm.normalize();
                gui.canvas.push_obj(&computations_cubic.V_norm,false);


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
            gui.canvas.pop(&computations_cubic.dual_m);
            gui.canvas.pop(&computations_cubic.m_grid);


            if(gui.InsideError.isChecked())
            {
                computations_cubic.f=get_scalar_field(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),0,0);
                computations_cubic.V=compute_field(computations_cubic.m,computations_cubic.f,gui.method.currentIndex());
                //computations_cubic.GT_grid=compute_values_on_grid(computations_cubic.m_grid, gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());

                bring_the_field_inside(computations_cubic.m,computations_cubic.m_grid,computations_cubic.V,computations_cubic.V_grid,gui.method.currentIndex());
                computations_cubic.GT=compute_ground_truth(computations_cubic.m_grid,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),LSDD);
                computations_cubic.err=estimate_error(computations_cubic.GT,computations_cubic.V_grid,computations_cubic.m_grid,gui.mode.currentIndex(),LSDD,gui.ErrType.currentIndex(),gui.type_of_vertices.currentIndex());
                computations_cubic.max=find_max_norm(computations_cubic.GT);
                computations_cubic.err_norm=heat_map_normalization(computations_cubic.err,0,computations_cubic.max,gui.sl_error_neg.value(),gui.sl_error_pos.value());
                computations_cubic.err_norm.copy_to_mesh(computations_cubic.m_grid);

                computations_cubic.m_grid.slice(gui.s);

                gui.canvas.pop(&computations_cubic.m);
                gui.canvas.push_obj(&computations_cubic.m_grid,false);
                computations_cubic.m_grid.show_in_wireframe(false);
                computations_cubic.m_grid.show_out_wireframe(false);
                computations_cubic.m_grid.show_in_texture1D(TEXTURE_1D_HSV);
                computations_cubic.m.show_in_vert_color();
                computations_cubic.m.show_out_vert_color();
                computations_cubic.m.slice(gui.s);
                gui.canvas.push_obj(&computations_cubic.m,false);
            }else{




                computations_cubic.f=get_scalar_field(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),0,0);
                computations_cubic.GT=compute_ground_truth(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
                computations_cubic.V=compute_field(computations_cubic.m,computations_cubic.f,gui.method.currentIndex());
                computations_cubic.err=estimate_error(computations_cubic.GT,computations_cubic.V,computations_cubic.m,gui.mode.currentIndex(),gui.method.currentIndex(),gui.ErrType.currentIndex(),gui.type_of_vertices.currentIndex());
                computations_cubic.max=find_max_norm(computations_cubic.GT);
                computations_cubic.err_norm=heat_map_normalization(computations_cubic.err,0,computations_cubic.max,gui.sl_error_neg.value(),gui.sl_error_pos.value());

                if(gui.method.currentIndex()==0)
                {
                    dual_mesh(computations_cubic.m, computations_cubic.dual_m,true);
                    computations_cubic.dual_m.slice(gui.s);
                    computations_cubic.err_norm.copy_to_mesh(computations_cubic.dual_m);

                    gui.canvas.pop(&computations_cubic.m);
                    gui.canvas.push_obj(&computations_cubic.dual_m,false);
                    computations_cubic.dual_m.show_in_wireframe(false);
                    computations_cubic.dual_m.show_out_wireframe(false);

                    computations_cubic.dual_m.show_in_texture1D(TEXTURE_1D_HSV);
                    computations_cubic.dual_m.show_out_texture1D(TEXTURE_1D_HSV);
                    gui.canvas.push_obj(&computations_cubic.m,false);
                }else
                {
                    computations_cubic.err_norm.copy_to_mesh(computations_cubic.m);
                    computations_cubic.m.show_in_texture1D(TEXTURE_1D_HSV);


                }
            }
            if(gui.but.isChecked())
            {
                computations_cubic.f=get_scalar_field(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),0,0);
                computations_cubic.V=compute_field(computations_cubic.m,computations_cubic.f,gui.method.currentIndex());
                computations_cubic.V_norm=computations_cubic.V;
                computations_cubic.V_norm.normalize();
                gui.canvas.push_obj(&computations_cubic.V_norm,false);


            }
            if(gui.ground_truth.isChecked())
            {
                computations_cubic.GT=compute_ground_truth(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
                computations_cubic.GT_norm=computations_cubic.GT;
                computations_cubic.GT_norm.normalize();
                gui.canvas.push_obj(&computations_cubic.GT_norm,false);

            }

        }
        gui.canvas.updateGL();


    });
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    QCheckBox::connect(&gui.InsideError, &QCheckBox::stateChanged, [&]()
    {
        gui.canvas.pop(&computations_cubic.dual_m);
        gui.canvas.pop(&computations_cubic.m_grid);
        //computations_cubic.dual_m.clear();


        if(gui.InsideError.isChecked())
        {


            computations_cubic.f=get_scalar_field(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),0,0);
            computations_cubic.V=compute_field(computations_cubic.m,computations_cubic.f,gui.method.currentIndex());


            bring_the_field_inside(computations_cubic.m,computations_cubic.m_grid,computations_cubic.V,computations_cubic.V_grid,gui.method.currentIndex());
            computations_cubic.GT=compute_ground_truth(computations_cubic.m_grid,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),LSDD);
            computations_cubic.err=estimate_error(computations_cubic.GT,computations_cubic.V_grid,computations_cubic.m_grid,gui.mode.currentIndex(),LSDD,gui.ErrType.currentIndex(),gui.type_of_vertices.currentIndex());
            computations_cubic.max=find_max_norm(computations_cubic.GT);
            computations_cubic.err_norm=heat_map_normalization(computations_cubic.err,0,computations_cubic.max,gui.sl_error_neg.value(),gui.sl_error_pos.value());
            computations_cubic.err_norm.copy_to_mesh(computations_cubic.m_grid);

            gui.canvas.pop(&computations_cubic.m);
            gui.canvas.push_obj(&computations_cubic.m_grid,false);
            computations_cubic.m_grid.slice(gui.s);
            computations_cubic.m_grid.show_in_wireframe(false);
            computations_cubic.m_grid.show_out_wireframe(false);
            computations_cubic.m_grid.show_in_texture1D(TEXTURE_1D_HSV);
            computations_cubic.m_grid.show_out_texture1D(TEXTURE_1D_HSV);
            computations_cubic.m.show_in_vert_color();
            computations_cubic.m.show_out_vert_color();

            gui.canvas.push_obj(&computations_cubic.m,false);



        }else{


            computations_cubic.f=get_scalar_field(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),0,0);
            computations_cubic.GT=compute_ground_truth(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
            computations_cubic.V=compute_field(computations_cubic.m,computations_cubic.f,gui.method.currentIndex());
            computations_cubic.err=estimate_error(computations_cubic.GT,computations_cubic.V,computations_cubic.m,gui.mode.currentIndex(),gui.method.currentIndex(),gui.ErrType.currentIndex(),gui.type_of_vertices.currentIndex());
            computations_cubic.max=find_max_norm(computations_cubic.GT);
            computations_cubic.err_norm=heat_map_normalization(computations_cubic.err,0,computations_cubic.max,gui.sl_error_neg.value(),gui.sl_error_pos.value());

            if(gui.method.currentIndex()==0)
            {
                computations_cubic.err_norm.copy_to_mesh(computations_cubic.dual_m);
                computations_cubic.dual_m.slice(gui.s);
                gui.canvas.pop(&computations_cubic.m);
                gui.canvas.push_obj(&computations_cubic.dual_m,false);
                computations_cubic.dual_m.show_in_wireframe(false);
                computations_cubic.dual_m.show_in_texture1D(TEXTURE_1D_HSV);
                gui.canvas.push_obj(&computations_cubic.m,false);



            }else
            {
                computations_cubic.err_norm.copy_to_mesh(computations_cubic.m);
                computations_cubic.m.show_in_texture1D(TEXTURE_1D_HSV);


            }





        }

        gui.canvas.updateGL();
    });
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    QComboBox::connect(&gui.type_of_vertices, static_cast<void(QComboBox::*)(int)>(&QComboBox::activated), [&]()
    {
        gui.canvas.pop(&computations_cubic.dual_m);
        gui.canvas.pop(&computations_cubic.m_grid);
        if(gui.InsideError.isChecked())
        {

            computations_cubic.f=get_scalar_field(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),0,0);
            computations_cubic.V=compute_field(computations_cubic.m,computations_cubic.f,gui.method.currentIndex());
            //computations_cubic.GT_grid=compute_values_on_grid(computations_cubic.m_grid, gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());

            bring_the_field_inside(computations_cubic.m,computations_cubic.m_grid,computations_cubic.V,computations_cubic.V_grid,gui.method.currentIndex());
            computations_cubic.GT=compute_ground_truth(computations_cubic.m_grid,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),LSDD);
            computations_cubic.err=estimate_error(computations_cubic.GT,computations_cubic.V_grid,computations_cubic.m_grid,gui.mode.currentIndex(),LSDD,gui.ErrType.currentIndex(),gui.type_of_vertices.currentIndex());
            computations_cubic.max=find_max_norm(computations_cubic.GT);
            computations_cubic.err_norm=heat_map_normalization(computations_cubic.err,0,computations_cubic.max,gui.sl_error_neg.value(),gui.sl_error_pos.value());
            computations_cubic.err_norm.copy_to_mesh(computations_cubic.m_grid);

            gui.canvas.pop(&computations_cubic.m);
            gui.canvas.push_obj(&computations_cubic.m_grid,false);
            computations_cubic.m_grid.show_in_wireframe(false);
            computations_cubic.m_grid.show_out_wireframe(false);
            computations_cubic.m_grid.show_in_texture1D(TEXTURE_1D_HSV);
            computations_cubic.m.show_in_vert_color();
            computations_cubic.m.show_out_vert_color();
            gui.canvas.push_obj(&computations_cubic.m,false);


        }else

        {
            computations_cubic.f=get_scalar_field(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),0,0);
            computations_cubic.GT=compute_ground_truth(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
            computations_cubic.V=compute_field(computations_cubic.m,computations_cubic.f,gui.method.currentIndex());
            computations_cubic.err=estimate_error(computations_cubic.GT,computations_cubic.V,computations_cubic.m,gui.mode.currentIndex(),gui.method.currentIndex(),gui.ErrType.currentIndex(),gui.type_of_vertices.currentIndex());
            computations_cubic.max=find_max_norm(computations_cubic.GT);
            computations_cubic.err_norm=heat_map_normalization(computations_cubic.err,0,computations_cubic.max,gui.sl_error_neg.value(),gui.sl_error_pos.value());

            if(gui.method.currentIndex()==0)
            {
                computations_cubic.err_norm.copy_to_mesh(computations_cubic.dual_m);
                gui.canvas.pop(&computations_cubic.m);
                gui.canvas.push_obj(&computations_cubic.dual_m,false);
                computations_cubic.dual_m.show_in_wireframe(false);
                computations_cubic.dual_m.show_in_texture1D(TEXTURE_1D_HSV);
                gui.canvas.push_obj(&computations_cubic.m,false);


            }else
            {


                computations_cubic.err_norm.copy_to_mesh(computations_cubic.m);
                computations_cubic.m.show_in_texture1D(TEXTURE_1D_HSV);


            }
        }
        gui.canvas.updateGL();

    });
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    QComboBox::connect(&gui.ErrType, static_cast<void(QComboBox::*)(int)>(&QComboBox::activated), [&]()
    {

        gui.canvas.pop(&computations_cubic.dual_m);
        gui.canvas.pop(&computations_cubic.m_grid);
        //computations_cubic.dual_m.clear();


        if(gui.InsideError.isChecked())
        {

            computations_cubic.f=get_scalar_field(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),0,0);
            computations_cubic.V=compute_field(computations_cubic.m,computations_cubic.f,gui.method.currentIndex());
            //computations_cubic.GT_grid=compute_values_on_grid(computations_cubic.m_grid, gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());

            bring_the_field_inside(computations_cubic.m,computations_cubic.m_grid,computations_cubic.V,computations_cubic.V_grid,gui.method.currentIndex());
            computations_cubic.GT=compute_ground_truth(computations_cubic.m_grid,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),LSDD);
            computations_cubic.err=estimate_error(computations_cubic.GT,computations_cubic.V_grid,computations_cubic.m_grid,gui.mode.currentIndex(),LSDD,gui.ErrType.currentIndex(),gui.type_of_vertices.currentIndex());
            computations_cubic.max=find_max_norm(computations_cubic.GT);
            computations_cubic.err_norm=heat_map_normalization(computations_cubic.err,0,computations_cubic.max,gui.sl_error_neg.value(),gui.sl_error_pos.value());
            computations_cubic.err_norm.copy_to_mesh(computations_cubic.m_grid);

            gui.canvas.pop(&computations_cubic.m);
            gui.canvas.push_obj(&computations_cubic.m_grid,false);
            computations_cubic.m_grid.show_in_wireframe(false);
            computations_cubic.m_grid.show_out_wireframe(false);
            computations_cubic.m_grid.show_in_texture1D(TEXTURE_1D_HSV);
            computations_cubic.m.show_in_vert_color();
            computations_cubic.m.show_out_vert_color();
            gui.canvas.push_obj(&computations_cubic.m,false);



        }else{
            computations_cubic.f=get_scalar_field(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),0,0);
            computations_cubic.GT=compute_ground_truth(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
            computations_cubic.V=compute_field(computations_cubic.m,computations_cubic.f,gui.method.currentIndex());
            computations_cubic.err=estimate_error(computations_cubic.GT,computations_cubic.V,computations_cubic.m,gui.mode.currentIndex(),gui.method.currentIndex(),gui.ErrType.currentIndex(),gui.type_of_vertices.currentIndex());
            computations_cubic.max=find_max_norm(computations_cubic.GT);
            computations_cubic.err_norm=heat_map_normalization(computations_cubic.err,0,computations_cubic.max,gui.sl_error_neg.value(),gui.sl_error_pos.value());

            if(gui.method.currentIndex()==0)
            {

                dual_mesh(computations_cubic.m, computations_cubic.dual_m,true);

                computations_cubic.err_norm.copy_to_mesh(computations_cubic.dual_m);

                gui.canvas.pop(&computations_cubic.m);
                gui.canvas.push_obj(&computations_cubic.dual_m,false);
                computations_cubic.dual_m.show_in_wireframe(false);
                computations_cubic.dual_m.show_out_wireframe(false);

                computations_cubic.dual_m.show_in_texture1D(TEXTURE_1D_HSV);
                gui.canvas.push_obj(&computations_cubic.m,false);
            }else
            {
                computations_cubic.err_norm.copy_to_mesh(computations_cubic.m);
                computations_cubic.m.show_in_texture1D(TEXTURE_1D_HSV);


            }

        }

        if(gui.but.isChecked())
        {
            computations_cubic.f=get_scalar_field(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),0,0);
            computations_cubic.V=compute_field(computations_cubic.m,computations_cubic.f,gui.method.currentIndex());
            computations_cubic.V_norm=computations_cubic.V;
            computations_cubic.V_norm.normalize();
            gui.canvas.push_obj(&computations_cubic.V_norm,false);


        }
        if(gui.ground_truth.isChecked())
        {
            computations_cubic.GT=compute_ground_truth(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
            computations_cubic.GT_norm=computations_cubic.GT;
            computations_cubic.GT_norm.normalize();
            gui.canvas.push_obj(&computations_cubic.GT_norm,false);

        }


        gui.canvas.updateGL();

    });
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    QComboBox::connect(&gui.mode, static_cast<void(QComboBox::*)(int)>(&QComboBox::activated), [&]()
    {


        gui.canvas.pop(&computations_cubic.dual_m);
        gui.canvas.pop(&computations_cubic.m_grid);
        //computations_cubic.dual_m.clear();


        if(gui.InsideError.isChecked())
        {
            computations_cubic.f=get_scalar_field(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),0,0);
            computations_cubic.V=compute_field(computations_cubic.m,computations_cubic.f,gui.method.currentIndex());
            //computations_cubic.GT_grid=compute_values_on_grid(computations_cubic.m_grid, gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex());

            bring_the_field_inside(computations_cubic.m,computations_cubic.m_grid,computations_cubic.V,computations_cubic.V_grid,gui.method.currentIndex());
            computations_cubic.GT=compute_ground_truth(computations_cubic.m_grid,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),LSDD);
            computations_cubic.err=estimate_error(computations_cubic.GT,computations_cubic.V_grid,computations_cubic.m_grid,gui.mode.currentIndex(),LSDD,gui.ErrType.currentIndex(),gui.type_of_vertices.currentIndex());
            computations_cubic.max=find_max_norm(computations_cubic.GT);
            computations_cubic.err_norm=heat_map_normalization(computations_cubic.err,0,computations_cubic.max,gui.sl_error_neg.value(),gui.sl_error_pos.value());
            computations_cubic.err_norm.copy_to_mesh(computations_cubic.m_grid);

            gui.canvas.pop(&computations_cubic.m);
            gui.canvas.push_obj(&computations_cubic.m_grid,false);
            computations_cubic.m_grid.show_in_wireframe(false);
            computations_cubic.m_grid.show_out_wireframe(false);
            computations_cubic.m_grid.show_in_texture1D(TEXTURE_1D_HSV);
            computations_cubic.m.show_in_vert_color();
            computations_cubic.m.show_out_vert_color();
            gui.canvas.push_obj(&computations_cubic.m,false);

        }else{


            computations_cubic.f=get_scalar_field(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),0,0);
            computations_cubic.GT=compute_ground_truth(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
            computations_cubic.V=compute_field(computations_cubic.m,computations_cubic.f,gui.method.currentIndex());
            computations_cubic.err=estimate_error(computations_cubic.GT,computations_cubic.V_grid,computations_cubic.m_grid,gui.mode.currentIndex(),LSDD,gui.ErrType.currentIndex(),gui.type_of_vertices.currentIndex());
            computations_cubic.max=find_max_norm(computations_cubic.GT);
            computations_cubic.err_norm=heat_map_normalization(computations_cubic.err,0,computations_cubic.max,gui.sl_error_neg.value(),gui.sl_error_pos.value());

            if(gui.method.currentIndex()==0)
            {
                computations_cubic.err_norm.copy_to_mesh(computations_cubic.dual_m);
                gui.canvas.pop(&computations_cubic.m);
                gui.canvas.push_obj(&computations_cubic.dual_m,false);
                computations_cubic.dual_m.show_in_wireframe(false);
                computations_cubic.dual_m.show_in_texture1D(TEXTURE_1D_HSV);
                gui.canvas.push_obj(&computations_cubic.m,false);


            }else
            {
                computations_cubic.err_norm.copy_to_mesh(computations_cubic.m);
                computations_cubic.m.show_in_texture1D(TEXTURE_1D_HSV);


            }
        }

        if(gui.but.isChecked())
        {
            computations_cubic.f=get_scalar_field(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),0,0);
            computations_cubic.V=compute_field(computations_cubic.m,computations_cubic.f,gui.method.currentIndex());
            computations_cubic.V_norm=computations_cubic.V;
            computations_cubic.V_norm.normalize();
            gui.canvas.push_obj(&computations_cubic.V_norm,false);


        }
        if(gui.ground_truth.isChecked())
        {
            computations_cubic.GT=compute_ground_truth(computations_cubic.m,gui.a.value(),gui.b.value(),gui.c.value(),gui.cb_function.currentIndex(),gui.method.currentIndex());
            computations_cubic.GT_norm=computations_cubic.GT;
            computations_cubic.GT_norm.normalize();
            gui.canvas.push_obj(&computations_cubic.GT_norm,false);

        }


        gui.canvas.updateGL();
    });

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    QPushButton::connect(&gui.load, &QPushButton::clicked, [&]()
    {
        std::string filename = QFileDialog::getOpenFileName(NULL, "Load mesh", ".", "3D Meshes (*.MESH *.HEDRA *.VTU *.VTK)").toStdString();
        if (!filename.empty())
        {
            DrawableTetmesh<> m(filename.c_str());
            computations_cubic.m=m;
            gui.canvas.push_obj( &computations_cubic.m);

        }
    });

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    QPushButton::connect(&gui.save, &QPushButton::clicked, [&]()
    {
//        ScalarField data=computations_cubic.err;
//        std::ofstream outfile;


//        outfile.open("/Users/Dala/Documents/MATLAB/Gradient_fields/" + std::to_string(10) + "result.csv");
//        for(int i=0; i<computations_cubic.m.num_verts();++i)
//        {
//            outfile<<data[i]<<","<<std::endl;
//        }

        std::string filename = QFileDialog::getSaveFileName(NULL, "Save mesh", ".", "3D Meshes (*.MESH *.TET *.HEDRA *.VTU *.VTK)").toStdString();
        if (!filename.empty()) computations_cubic.m.save(filename.c_str());


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

        computations_cubic.dual_m.show_in_vert_color();
        gui.canvas.pop_all_occurrences_of(DRAWABLE_VECTOR_FIELD);
        gui.canvas.pop(&computations_cubic.dual_m);
        gui.canvas.updateGL();

    });
}






