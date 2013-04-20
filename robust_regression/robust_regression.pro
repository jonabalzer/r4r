TEMPLATE = app
CONFIG += console
CONFIG -= qt
TARGET = lmtest

SOURCES += main.cpp

HEADERS += \
    main.h

QMAKE_CXXFLAGS += -std=c++0x

LIBS += -L$$OUT_PWD/../r4r_core/ -lr4r_core

INCLUDEPATH += $$PWD/../r4r_core
DEPENDPATH += $$PWD/../r4r_core

