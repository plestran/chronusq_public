set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I /usr/include/python2.6")

add_library(pythonapi SHARED singleslater_api.cpp molecule_api.cpp 
aointegrals_api.cpp chronusq_api.cpp classtools_api.cpp fileio_api.cpp
tools_api.cpp basisset_api.cpp)

target_link_libraries(pythonapi LINK_PUBLIC basisset classtools controls
mointegrals molecule tools sdresponse realtime mollerplesset singleslater 
aointegrals worker quasinewton grid fileio cerr gauinterface ${LIBINT_LIB} m 
boost_python python2.6 lapack blas)

