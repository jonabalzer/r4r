TEMPLATE = app
CONFIG += console
TARGET = lmtest

SOURCES += main.cpp

HEADERS += main.h

QMAKE_CXXFLAGS += -std=c++0x

DEPENDPATH += $$PWD/../r4r_core

INCLUDEPATH += $$PWD/../r4r_core

unix:!symbian|win32 {

    # include paths
    INCLUDEPATH += /usr/include/r4r/

    LIBS += -L$$OUT_PWD/../r4r_core/ \
            -lr4r_core

    target.path = $$OUT_PWD/../bin
    INSTALLS += target \

}
