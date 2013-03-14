TEMPLATE = app
CONFIG += console
CONFIG -= qt
QMAKE_CXXFLAGS += -std=c++0x

SOURCES += main.cpp

DEFINES += HAVE_FFTW
unix:!macx:!symbian: LIBS += -L$$OUT_PWD/../r4r_core/ -lr4r_core

INCLUDEPATH += $$PWD/../r4r_core
DEPENDPATH += $$PWD/../r4r_core

unix:!macx:!symbian: LIBS += -L$$OUT_PWD/../r4r_motion/ -lr4r_motion

INCLUDEPATH += $$PWD/../r4r_motion
DEPENDPATH += $$PWD/../r4r_motion

LIBS += -L/usr/local/lib \
     -lopencv_core\
     -lopencv_highgui\
     -lopencv_video\
     -lopencv_imgproc\
     -lopencv_features2d\
     -lopencv_calib3d
