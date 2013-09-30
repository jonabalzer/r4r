#-------------------------------------------------
#
# Project created by QtCreator 2013-09-29T19:35:00
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = tvdenoising
TEMPLATE = app

QMAKE_CXXFLAGS += -std=c++0x -O3

SOURCES += main.cpp\
           mainwindow.cpp \
           nabla.cpp

HEADERS  += mainwindow.h \
            nabla.h

FORMS    += mainwindow.ui

INCLUDEPATH += $$PWD/../r4r_core

DEPENDPATH += $$PWD/../r4r_core

target.path = $$OUT_PWD/../bin
INSTALLS += target

unix:!symbian|win32: {

    LIBS += -L$$OUT_PWD/../r4r_core/ \
            -lr4r_core \

   INCLUDEPATH += /usr/include/r4r/

}

RESOURCES += \
    icon.qrc
