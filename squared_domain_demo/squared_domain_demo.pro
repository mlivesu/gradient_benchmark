LIBS_PATH       = #Enter your path to cinolib folder here
TEMPLATE        = app
TARGET          = squared_domain_demo
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

## enable Triangle
DEFINES     += CINOLIB_USES_TRIANGLE
INCLUDEPATH += $$LIBS_PATH
LIBS        += -L$$LIBS_PATH -ltriangle



SOURCES += main.cpp \
    triangulations.cpp
SOURCES += computations.cpp
SOURCES += picking.cpp
SOURCES += gui.cpp

HEADERS += gui.h \
    triangulations.h
HEADERS += computations.h
HEADERS += functions.h
HEADERS += picking.h
