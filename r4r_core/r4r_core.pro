#-------------------------------------------------
#
# Project created by QtCreator 2013-02-21T16:20:49
#
#-------------------------------------------------

QMAKE_CXXFLAGS += -std=c++0x -fopenmp -msse2

TARGET = r4r_core
TEMPLATE = lib

DEFINES += R4R_CORE_LIBRARY

SOURCES += \
    trafo.cpp \
    sarray.cpp \
    rect.cpp \
    precond.cpp \
    params.cpp \
    lm.cpp \
    kfilter.cpp \
    iter.cpp \
    intimg.cpp \
    interp.cpp \
    factor.cpp \
    darray.cpp \
    cam.cpp \
    pegasos.cpp \
    rutils.cpp

HEADERS += \
    types.h \
    trafo.h \
    sarray.h \
    rect.h \
    precond.h \
    params.h \
    lm.h \
    kfilter.h \
    iter.h \
    intimg.h \
    interp.h \
    factor.h \
    darray.h \
    cam.h \
    pegasos.h \
    rutils.h

LIBS += -L/usr/local/lib \
     -lopencv_core \
     -lopencv_highgui \
     -lopencv_video \
     -lopencv_imgproc \
     -lopencv_features2d \
     -lopencv_calib3d \
     -llapack \
     -lgomp

FORTRAN_SOURCES += \
    spdot.f \
    intkernel.f \
    idkernel.f \
    helkernel.f \
    chikernel.f

fortran.output = ${QMAKE_FILE_BASE}.o
fortran.commands = f95 -c ${QMAKE_FILE_NAME} -o ${QMAKE_FILE_OUT} -fPIC
fortran.input = FORTRAN_SOURCES
QMAKE_EXTRA_COMPILERS += fortran
