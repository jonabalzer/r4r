#-------------------------------------------------
#
# Project created by QtCreator 2013-02-21T19:42:29
#
#-------------------------------------------------

QT       += core gui

TARGET = mv_descriptor_learning
TEMPLATE = app


QMAKE_CXXFLAGS += -std=c++0x -O3

SOURCES += main.cpp\
           mainwindow.cpp \
           preferences.cpp

HEADERS  += mainwindow.h \
            preferences.h

FORMS    += mainwindow.ui \
            preferences.ui

unix:!symbian|win32: LIBS += -L/usr/local/lib/\
                             -lopencv_core\
                             -lopencv_highgui\
                             -lopencv_video\
                             -lopencv_imgproc\
                             -lopencv_features2d\
                             -lopencv_calib3d \
                             -llapack \
                             -lgomp


win32:CONFIG(release, debug|release): LIBS += -L$$OUT_PWD/../r4r_core/release/ -lr4r_core
else:win32:CONFIG(debug, debug|release): LIBS += -L$$OUT_PWD/../r4r_core/debug/ -lr4r_core
else:unix:!symbian: LIBS += -L$$OUT_PWD/../r4r_core/ -lr4r_core

INCLUDEPATH += $$PWD/../r4r_core
DEPENDPATH += $$PWD/../r4r_core

win32:CONFIG(release, debug|release): LIBS += -L$$OUT_PWD/../r4r_motion/release/ -lr4r_motion
else:win32:CONFIG(debug, debug|release): LIBS += -L$$OUT_PWD/../r4r_motion/debug/ -lr4r_motion
else:unix:!symbian: LIBS += -L$$OUT_PWD/../r4r_motion/ -lr4r_motion

INCLUDEPATH += $$PWD/../r4r_motion
DEPENDPATH += $$PWD/../r4r_motion
