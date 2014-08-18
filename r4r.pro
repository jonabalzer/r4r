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

TEMPLATE = subdirs

CONFIG += ordered

SUBDIRS += r4r_core \
           r4r_motion \
           r4r_hardware \
           r4r_reconstruction

# optionally build unit tests
equals( HAVE_UNITTESTS, 1 ) {

    SUBDIRS += unit_tests

}

# optionally include examples in the build
equals( HAVE_EXAMPLES, 1 ) {

    SUBDIRS += robust_regression \
               tvdenoising \
               mv_descriptor_learning \
               clam

}

QMAKE_CXXFLAGS += -std=c++0x

# doxygen target
dox.target = doc
dox.commands = doxygen $$PWD/doxyfile
QMAKE_EXTRA_TARGETS += dox
