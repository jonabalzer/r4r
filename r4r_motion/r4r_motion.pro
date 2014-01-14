######################################################################################
#
# Copyright (c) 2013, Jonathan Balzer
#
# All rights reserved.
#
# This file is part of the R4R library.
#
# The R4R library is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# The R4R library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with the R4R library. If not, see <http://www.gnu.org/licenses/>.
#
######################################################################################

QMAKE_CXXFLAGS += -std=c++0x -O3

TARGET = r4r_motion
TEMPLATE = lib

DEFINES += R4R_MOTION_LIBRARY

SOURCES += tracker.cpp \
    stracker.cpp \
    mtracker.cpp \
    lk.cpp \
    feature.cpp \
    descriptor.cpp \
    basic.cpp \
    tracklet.cpp \
    dagg.cpp \
    descspecial.cpp

HEADERS += \
    tracker.h \
    stracker.h \
    mtracker.h \
    lk.h \
    feature.h \
    descriptor.h \
    basic.h \
    tracklet.h \
    dagg.h \
    descspecial.h

# add last flag depending on whether FFTW is present
packagesExist(fftw3) {
    DEFINES += HAVE_FFTW
}

# make sure that r4r_core is up to date
DEPENDPATH += $$PWD/../r4r_core

# local inlude path
INCLUDEPATH += $$PWD/../r4r_core

# since we have dependencies, create a pc file
CONFIG += create_prl no_install_prl create_pc

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
            -lopencv_video \
            -L$$OUT_PWD/../r4r_core \
            -lr4r_core

    contains(DEFINES,HAVE_FFTW) {

        LIBS += -lfftw3

    }

}


