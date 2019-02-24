LIBS_PATH       = $$PWD/../../lib
TEMPLATE        = app
TARGET          = batch_program_cubic_domain
QT             += core opengl
CONFIG         += c++11
CONFIG         -= app_bundle
INCLUDEPATH    += $$LIBS_PATH/CinoLib/include
INCLUDEPATH    += $$LIBS_PATH/CinoLib/external/eigen
DEFINES        += CINOLIB_USES_OPENGL
DEFINES        += CINOLIB_USES_QT
QMAKE_CXXFLAGS += -Wno-deprecated-declarations # gluQuadric gluSphere and gluCylinde are deprecated in macOS 10.9

# just for Linux
unix:!macx {
DEFINES += GL_GLEXT_PROTOTYPES
LIBS    += -lGLU
}

## enable Tetgen
DEFINES     += CINOLIB_USES_TETGEN
DEFINES     += TETLIBRARY
INCLUDEPATH += $$LIBS_PATH
LIBS        += -L$$LIBS_PATH -ltet

## enable Triangle
DEFINES     += CINOLIB_USES_TRIANGLE
INCLUDEPATH += $$LIBS_PATH
LIBS        += -L$$LIBS_PATH -ltriangle



SOURCES += \
        main.cpp \
    computations_cubic.cpp \
    polygonsoup.cpp \
    meshes.cpp

HEADERS += \
    computations_cubic.h \
    functions_cubic.h \
    polygonsoup.h \
    meshes.h


