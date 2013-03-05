TEMPLATE = subdirs

SUBDIRS += \
    r4r_core \
    r4r_motion \
    robust_regression \
    mv_descriptor_learning \
    scratch



QMAKE_CXXFLAGS += -std=c++0x
