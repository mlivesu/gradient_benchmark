#ifndef GUI_H
#define GUI_H

#include <QWidget>
#include <QComboBox>
#include <QGroupBox>
#include <QPushButton>
#include<QRadioButton>
#include <QSpinBox>
#include<QCheckBox>
#include <QSlider>
#include <cinolib/gui/qt/glcanvas.h>
#include <cinolib/meshes/meshes.h>
#include <cinolib/meshes/mesh_slicer.h>
#include <cinolib/scalar_field.h>
#include <QFileDialog>
#include <cinolib/drawable_vector_field.h>
#include <cinolib/gradient.h>
#include "polygonsoup.h"
#include "computations_cubic.h"
#include "meshes.h"





using namespace cinolib;
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

typedef struct
{
    QWidget           widget;
    cinolib::GLcanvas canvas;
    QGroupBox         error;
    QGroupBox         choose_heat_map;
    QGroupBox         show;
    QGroupBox         slicing;
    QComboBox         cb_function;
    QComboBox         method;
    QComboBox         choose_tri;
    QCheckBox         but;
    QCheckBox         ground_truth;
    QComboBox         type_of_vertices;
    QCheckBox         wireframe;
    QCheckBox         show_heat_map;
    QCheckBox         cb_slice_flip_x;
    QCheckBox         cb_slice_flip_y;
    QCheckBox         cb_slice_flip_z;
    QPushButton       reset;
    QRadioButton      scalar_field;
    QRadioButton      Err;
    QCheckBox         InsideError;
    QPushButton       save;
    QPushButton       load;
    QSpinBox          N;
    QSpinBox          anis;
    QSpinBox          a;
    QSpinBox          b;
    QSpinBox          c;
    QComboBox         ErrType;
    QComboBox         mode;
    QSlider           sl_scalar_field;
    QSlider           sl_error;
    QSlider           sl_vector_fields;
    QSlider           sl_slice_x;
    QSlider           sl_slice_y;
    QSlider           sl_slice_z;
    SlicerState       s;

}
GUI;

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void init_gui   (GUI & gui);
void init_events(GUI & gui);

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#endif // GUI_H
