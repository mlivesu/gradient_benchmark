#include <QApplication>
#include <cinolib/gui/qt/qt_gui_tools.h>
#include <cinolib/point_inside_mesh.h>
#include <iostream>
#include "gui.h"

using namespace cinolib;

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

COMPUTATIONS computations;

int main(int argc, char **argv)
{
    QApplication a(argc, argv);

    GUI gui;
    init_gui(gui);
    init_events(gui);

    make_triangulation(computations.m,1,25);
    make_grid(computations.m_grid,100);
    gui.canvas.push_obj(&computations.m);

    // CMD+1 to show mesh controls.
    SurfaceMeshControlPanel<DrawableTrimesh<>> panel(&computations.m, &gui.canvas);
    QApplication::connect(new QShortcut(QKeySequence(Qt::CTRL+Qt::Key_1), &gui.canvas), &QShortcut::activated, [&](){panel.show();});

    return a.exec();
}
