
set(CMAKE_CXX_FLAGS "-Wall -std=c++0x ${CMAKE_CXX_FLAGS} -fopenmp -O3") 

set(SOURCES_H
    cam.h
    darray.h
    factor.h
    interp.h
    intimg.h
    iter.h
    kernels.h
    kfilter.h
    lm.h
    params.h
    pegasos.h
    precond.h
    rect.h
    sarray.h
    trafo.h
    types.h
    rutils.h
    vecn.h)
    
set(SOURCES_CPP
    cam.cpp
    darray.cpp
    factor.cpp
    interp.cpp
    intimg.cpp
    iter.cpp
    kernels.cpp
    kfilter.cpp
    lm.cpp
    params.cpp
    pegasos.cpp
    precond.cpp
    rect.cpp
    sarray.cpp
    trafo.cpp
    rutils.cpp
    vecn.cpp)
    
add_library(r4r_core SHARED ${SOURCES_CPP})
install(DIRECTORY DESTINATION include/r4r/r4r_core)
install(FILES ${SOURCES_H} DESTINATION include/r4r/r4r_core/)
install(TARGETS r4r_core LIBRARY DESTINATION lib)
