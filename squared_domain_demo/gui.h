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
#include <cinolib/scalar_field.h>
#include <cinolib/drawable_vector_field.h>
#include <QFileDialog>
#include <GG_gradient.h>
#include "triangulations.h"
#include "computations.h"



using namespace cinolib;
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

typedef struct
{
    QWidget           widget;
    cinolib::GLcanvas canvas;
    QGroupBox         error;
    QGroupBox         choose_heat_map;
    QComboBox         cb_function;
    QComboBox         method;
    QComboBox         choose_tri;
    QCheckBox         but;
    QCheckBox         ground_truth;
    QComboBox         type_of_vertices;
    QCheckBox         wireframe;
    QCheckBox         show_heat_map;
    QPushButton       reset;
    QRadioButton      scalar_field;
    QRadioButton      Err;
    QCheckBox         InsideError;
    QCheckBox         boundaries;
    QPushButton       save;
    QSpinBox          N;
    QSpinBox          anis;
    QSpinBox          a;
    QSpinBox          b;
    QSpinBox          c;
    QComboBox         ErrType;
    QComboBox         mode;
    QSlider           sl_scalar_field;
    QSlider           sl_error_neg;
    QSlider           sl_error_pos;
    QSlider           sl_vector_fields;
}
GUI;

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void init_gui   (GUI & gui);
void init_events(GUI & gui);

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#endif // GUI_H
