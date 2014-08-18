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
    rbuffer.h \
    unionfind.h

unix:!symbian|win32 {

    # create pkg config files
    CONFIG += warn_off create_prl no_install_prl create_pc

    # find LAPACK and BLAS
    LAPACK = $$system(find /usr -name liblapack* 2>/dev/null)
    isEmpty(LAPACK) {
        error("Could not resolve dependency on LAPACK.")
    }
    else {
        LIBS += -llapack
        DEFINES += HAVE_LAPACK
    }
    BLAS = $$system(find /usr -name libblas* 2>/dev/null)
    isEmpty(BLAS) {
        warning("Could not resolve dependency on BLAS.")
    }
    else {

        LIBS += -lblas
        DEFINES += HAVE_BLAS

    }

    # find OpenCV
    packagesExist(opencv) {

        LIBS += -lopencv_core \
                -lopencv_imgproc

        # this is not good, because it will spawn too many dependencies
        #CONFIG += link_pkgconfig
        #PKGCONFIG += opencv

    } else {
       error("Could not resolve mandatory dependence on OpenCV...")

    }

    # find TBB
    #packagesExist(tbb) { DEFINES += HAVE_TBB }

    # find OpenEXR
#    packagesExist(OpenEXR) {

#        CONFIG += link_pkgconfig
#        PKGCONFIG += OpenEXR
#        DEFINES += HAVE_EXR

#    } else {
#        warning("Could not resolve dependence on OpenEXR...")
#    }

    # install target, FIXME: add target for pkg file
    headers.files = $$HEADERS
    headers.path = /usr/include/r4r/
    target.path = /usr/lib/

    INSTALLS += target \
                headers

}
