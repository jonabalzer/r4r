TEMPLATE = app
CONFIG += console
CONFIG -= qt
TARGET = lmtest

SOURCES += main.cpp

HEADERS += \
    main.h



win32:CONFIG(release, debug|release): LIBS += -L$$OUT_PWD/../r4r_core/release/ -lr4r_core
else:win32:CONFIG(debug, debug|release): LIBS += -L$$OUT_PWD/../r4r_core/debug/ -lr4r_core
else:unix:!symbian: LIBS += -L$$OUT_PWD/../r4r_core/ -lr4r_core

INCLUDEPATH += $$PWD/../r4r_core
DEPENDPATH += $$PWD/../r4r_core
