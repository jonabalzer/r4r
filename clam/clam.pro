#-------------------------------------------------
#
# Project created by QtCreator 2013-07-05T14:55:04
#
#-------------------------------------------------

QT       += core gui opengl

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = clam
TEMPLATE = app

QMAKE_CXXFLAGS += -std=c++0x -O3

SOURCES += main.cpp\
        mainwindow.cpp \
        preferences.cpp \
        viewer.cpp

HEADERS += mainwindow.h \
           preferences.h \
           viewer.h

FORMS += mainwindow.ui \
         preferences.ui

INCLUDEPATH += $$PWD/../r4r_core \
               $$PWD/../r4r_motion

DEPENDPATH += $$PWD/../r4r_core \
              $$PWD/../r4r_motion

target.path = $$OUT_PWD/../bin
INSTALLS += target

RESOURCES += clamicon.qrc

unix:!symbian|win32: {

    LIBS += -L/usr/local/lib/\
            -lopencv_core\
            -lopencv_highgui\
            -lopencv_video\
            -lopencv_imgproc\
            -lopencv_features2d\
            -lopencv_calib3d \
            -llapack \
            -lgomp \
            -L$$OUT_PWD/../r4r_core/ \
            -lr4r_core \
            -L$$OUT_PWD/../r4r_motion/ \
            -lr4r_motion

   INCLUDEPATH += /usr/include/r4r/

}

