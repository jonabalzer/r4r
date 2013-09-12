#-------------------------------------------------
#
# Project created by QtCreator 2013-02-21T19:42:29
#
#-------------------------------------------------

QT += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = mvdl
TEMPLATE = app

QMAKE_CXXFLAGS += -std=c++0x -O3

SOURCES += main.cpp\
           mainwindow.cpp \
           preferences.cpp

HEADERS += mainwindow.h \
           preferences.h

FORMS += mainwindow.ui \
         preferences.ui

INCLUDEPATH += $$PWD/../r4r_core \
               $$PWD/../r4r_motion

DEPENDPATH += $$PWD/../r4r_core \
              $$PWD/../r4r_motion

target.path = $$OUT_PWD/../bin
INSTALLS += target

unix:!symbian|win32: {

    LIBS += -L$$OUT_PWD/../r4r_core/ \
            -L$$OUT_PWD/../r4r_motion/ \
            -lr4r_core \
            -lr4r_motion \
            -L/usr/local/lib/\
            -lopencv_core\
            -lopencv_highgui\
            -lopencv_video\
            -lopencv_imgproc\
            -lopencv_features2d\
            -lopencv_calib3d \
            -llapack \
            -lgomp

   INCLUDEPATH += /usr/include/r4r/

}

