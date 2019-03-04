LIBS_PATH       = $$PWD/../../gradient_fields/lib
DATA_PATH       = $$PWD/../../data

TEMPLATE        = app
TARGET          = batch_program_squared_domain
QT             += core opengl
CONFIG         += c++11
CONFIG         -= app_bundle
INCLUDEPATH    += $$LIBS_PATH/CinoLib/include
INCLUDEPATH    += $$LIBS_PATH/CinoLib/external/eigen
DEFINES        += CINOLIB_USES_OPENGL
DEFINES        += CINOLIB_USES_QT
DEFINES        += DATA_PATH=$$DATA_PATH
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



SOURCES += \
        main.cpp \
    computations.cpp \
    triangulations.cpp

HEADERS += \
    computations.h \
    triangulations.h \
    functions.h
