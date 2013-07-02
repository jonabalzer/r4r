
set(CMAKE_CXX_FLAGS "-Wall -std=c++0x ${CMAKE_CXX_FLAGS}") 

set(SOURCES_H
    cam.h
    darray.h
    factor.h
    interp.h
    intimg.h
    iter.h
    kfilter.h
    lm.h
    params.h
    pegasos.h
    precond.h
    rect.h
    sarray.h
    trafo.h
    types.h
    rutils.h)
    
set(SOURCES_CPP
    cam.cpp
    darray.cpp
    factor.cpp
    interp.cpp
    intimg.cpp
    iter.cpp
    kfilter.cpp
    lm.cpp
    params.cpp
    pegasos.cpp
    precond.cpp
    rect.cpp
    sarray.cpp
    trafo.cpp
    rutils.cpp)
    
add_library(r4r_core SHARED ${SOURCES_CPP})
