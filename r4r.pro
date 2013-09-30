TEMPLATE = subdirs

CONFIG += ordered

SUBDIRS += r4r_core \
           r4r_motion \
           r4r_hardware \
    tvdenoising

equals( HAVE_EXAMPLES, 1 ) {

    SUBDIRS += robust_regression \
               mv_descriptor_learning \
               clam

}

QMAKE_CXXFLAGS += -std=c++0x

# doxygen target
dox.target = doc
dox.commands = doxygen $$PWD/doxyfile
QMAKE_EXTRA_TARGETS += dox
