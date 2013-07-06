
include(FindPkgConfig)
pkg_check_modules(FFTW3 REQUIRED fftw3)
add_definitions(-DHAVE_FFTW)
message(STATUS "FFTW3_LIBRARIES: ${FFTW3_LIBRARIES}")

set(CMAKE_CXX_FLAGS "-Wall -std=c++0x ${CMAKE_CXX_FLAGS} -fopenmp -O3") 

set(SOURCES_Hbasic.h
    dagg.h
    descriptor.h
    descspecial.h
    feature.h
    lk.h
    map.h
    mtracker.h
    stracker.h
    tracker.h
    tracklet.h)

set(SOURCES_CPP
    basic.cpp
    dagg.cpp
    descriptor.cpp
    descspecial.cpp
    feature.cpp
    lk.cpp
    map.cpp
    mtracker.cpp
    stracker.cpp
    tracker.cpp
    tracklet.cpp)
    
add_library(r4r_motion SHARED ${SOURCES_CPP})
target_link_libraries(r4r_motion ${FFTW3_LIBRARIES})
