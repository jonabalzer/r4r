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

QT += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = mvdl
TEMPLATE = app

QMAKE_CXXFLAGS += -std=c++0x -O3 -fopenmp

SOURCES += main.cpp\
           mainwindow.cpp \
           preferences.cpp

HEADERS += mainwindow.h \
           preferences.h

FORMS += mainwindow.ui \
         preferences.ui

INCLUDEPATH += $$PWD/../r4r_core \
               $$PWD/../r4r_motion

DEPENDPATH += $$PWD/../r4r_core \
              $$PWD/../r4r_motion

target.path = $$OUT_PWD/../bin
INSTALLS += target

unix:!symbian|win32: {

    LIBS += -L$$OUT_PWD/../r4r_core/ \
            -L$$OUT_PWD/../r4r_motion/ \
            -lr4r_core \
            -lr4r_motion \
            -L/usr/local/lib/\
            -lopencv_core\
            -lopencv_highgui\
            -lopencv_video\
            -lopencv_imgproc\
            -lopencv_features2d\
            -lopencv_calib3d \
            -llapack \
            -lgomp

   INCLUDEPATH += /usr/include/r4r/

}

