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

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = tvdenoising
TEMPLATE = app

QMAKE_CXXFLAGS += -std=c++0x -O2

SOURCES += main.cpp\
           mainwindow.cpp \
           nabla.cpp

HEADERS  += mainwindow.h \
            nabla.h

FORMS    += mainwindow.ui

INCLUDEPATH += $$PWD/../r4r_core

DEPENDPATH += $$PWD/../r4r_core

target.path = $$OUT_PWD/../bin
INSTALLS += target

packagesExist(tbb) {
    DEFINES += HAVE_TBB
    CONFIG(debug,debug|release):DEFINES += TBB_USE_DEBUG
}


unix:!symbian|win32: {

    LIBS += -L$$OUT_PWD/../r4r_core/ \
            -lr4r_core \

    contains(DEFINES,HAVE_TBB) {
            LIBS += -ltbb
            #CONFIG(release,debug|release):LIBS+= -ltbb
            #CONFIG(debug,debug|release):LIBS+= -ltbb_debug
    }

   INCLUDEPATH += /usr/include/r4r/

}

RESOURCES += \
    icon.qrc
