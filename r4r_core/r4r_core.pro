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

#QT -= core gui

QMAKE_CXXFLAGS += -std=c++0x -O3 -msse4

TARGET = r4r_core
TEMPLATE = lib

DEFINES += R4R_CORE_LIBRARY

SOURCES += \
    trafo.cpp \
    sarray.cpp \
    rect.cpp \
    precond.cpp \
    params.cpp \
    lm.cpp \
    kfilter.cpp \
    iter.cpp \
    intimg.cpp \
    interp.cpp \
    factor.cpp \
    darray.cpp \
    cam.cpp \
    pegasos.cpp \
    rutils.cpp \
    kernels.cpp \
    splinecurve.cpp \
    vecn.cpp \
    image.cpp \
    types.cpp

HEADERS += \
    types.h \
    trafo.h \
    sarray.h \
    rect.h \
    precond.h \
    params.h \
    lm.h \
    kfilter.h \
    iter.h \
    intimg.h \
    interp.h \
    factor.h \
    darray.h \
    cam.h \
    pegasos.h \
    rutils.h \
    kernels.h \
    splinecurve.h \
    vecn.h \
    image.h \
    rbuffer.h

# see if intel tbb existst
packagesExist(tbb) {
    DEFINES += HAVE_TBB
    CONFIG(debug,debug|release):DEFINES += TBB_USE_DEBUG
}

unix:!symbian|win32 {

    headers.files = $$HEADERS
    headers.path = /usr/include/r4r/

    target.path = /usr/lib/

    INSTALLS += target \
                headers

    # what about clean target?

    LIBS += -L/usr/local/lib \
            -lopencv_core \
            -lopencv_highgui \
            -lopencv_video \
            -lopencv_imgproc \
            -lopencv_features2d \
            -lopencv_calib3d \
            -llapack

    contains(DEFINES,HAVE_TBB) {
            LIBS += -ltbb
            #CONFIG(release,debug|release):LIBS+= -ltbb
            #CONFIG(debug,debug|release):LIBS+= -ltbb_debug
    }




}
