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

QT       += core gui opengl

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = clam
TEMPLATE = app

QMAKE_CXXFLAGS += -std=c++0x -O3

SOURCES += main.cpp\
        mainwindow.cpp \
        preferences.cpp

HEADERS += mainwindow.h \
           preferences.h

FORMS += mainwindow.ui \
         preferences.ui

INCLUDEPATH += $$PWD/../r4r_core \
               $$PWD/../r4r_motion \
               $$PWD/../r4r_reconstruction

DEPENDPATH += $$PWD/../r4r_core \
              $$PWD/../r4r_motion \
              $$PWD/../r4r_reconstruction

target.path = $$OUT_PWD/../bin
INSTALLS += target

RESOURCES += clamicon.qrc

# find OpenMesh library
OM = $$system(find /usr -name libOpenMeshCore* 2>/dev/null)
isEmpty(OM) {
    error("Could not resolve dependency on OpenMesh.")
}
OM = $$first(OM)
OMLIBPATH = $$dirname(OM)

unix:!symbian|win32: {

    LIBS += -L$$OUT_PWD/../r4r_core/ \
            -L$$OUT_PWD/../r4r_motion/ \
            -L$$OUT_PWD/../r4r_reconstruction/ \
            -lr4r_core \
            -lr4r_motion \
            -lr4r_reconstruction \
            -L/usr/local/lib/\
            -lopencv_core\
            -lopencv_highgui\
            -lopencv_video\
            -lopencv_imgproc\
            -lopencv_features2d\
            -lopencv_calib3d \
            -llapack \
            -lgomp \
            -L$$OMLIBPATH \
            -lOpenMeshCore \
            -lOpenMeshTools


   INCLUDEPATH += /usr/include/r4r/

}

