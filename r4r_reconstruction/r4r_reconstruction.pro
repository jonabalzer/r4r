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

QT += opengl

QMAKE_CXXFLAGS += -std=c++0x

TARGET = r4r_reconstruction
TEMPLATE = lib

DEFINES += R4R_RECONSTRUCTION_LIBRARY

SOURCES += viewer.cpp \
    trimesh.cpp \
    bbox.cpp \
    pcl.cpp

HEADERS += viewer.h \
    trimesh.h \
    bbox.h \
    pcl.h

# find OpenMesh library
OM = $$system(find /usr -name libOpenMeshCore* 2>/dev/null)
isEmpty(OM) {
    error("Could not resolve dependency on OpenMesh.")
}
OM = $$first(OM)
OMLIBPATH = $$dirname(OM)

LIBS += -L$$OMLIBPATH \
        -lOpenMeshCore \
        -lOpenMeshTools

# local inlude path
INCLUDEPATH += $$PWD/../r4r_core \
               $$PWD/../r4r_motion

# make sure that r4r_core is up to date
DEPENDPATH += $$PWD/../r4r_core \
              $$PWD/../r4r_motion

CONFIG += create_prl no_install_prl create_pc

unix:!symbian|win32 {

    headers.files = $$HEADERS
    headers.path = /usr/include/r4r/

    target.path = /usr/lib/

    INSTALLS += target \
                headers

}
