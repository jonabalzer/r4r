QMAKE_CXXFLAGS += -std=c++0x -O3

# add last flag depending on whether FFTW is present
DEFINES += HAVE_FFTW

TARGET = r4r_motion
TEMPLATE = lib

DEFINES += R4R_MOTION_LIBRARY

SOURCES += tracker.cpp \
    stracker.cpp \
    #ptracker.cpp \
    mtracker.cpp \
    lk.cpp \
    feature.cpp \
    descriptor.cpp \
    basic.cpp \
    #tsttrack.cpp \
    tracklet.cpp \
    dagg.cpp \
    descspecial.cpp

HEADERS += \
    tracker.h \
    stracker.h \
    #ptracker.h \
    mtracker.h \
    lk.h \
    feature.h \
    descriptor.h \
    basic.h \
    #tsttrack.h \
    tracklet.h \
    dagg.h \
    descspecial.h


# make sure that r4r_core is up to date
DEPENDPATH += $$PWD/../r4r_core

# local inlude path
INCLUDEPATH += $$PWD/../r4r_core

unix:!symbian|win32 {

    # set install path
    headers.files = $$HEADERS
    headers.path = /usr/include/r4r/
    target.path = /usr/lib/

    INSTALLS += target \
                headers

    # include paths
    INCLUDEPATH += /usr/include/r4r/

    LIBS += -L/usr/local/lib/ \
            -lopencv_features2d \
            -lopencv_video

    contains( DEFINES, HAVE_FFTW ) {

        LIBS += -lfftw3

    }

}


