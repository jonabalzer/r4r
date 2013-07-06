TEMPLATE = subdirs

SUBDIRS += \
    r4r_core \
    r4r_motion \
    robust_regression \
    mv_descriptor_learning \
    scratch \
    r4r_hardware \
    iotest

QMAKE_CXXFLAGS += -std=c++0x

# doxygen target
dox.target = doc
dox.commands = doxygen ../src/doxyfile
QMAKE_EXTRA_TARGETS += dox
