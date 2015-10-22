set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I ${PYTHON_INCLUDE_DIRS}")


add_library(pythonapi SHARED singleslater_api.cpp molecule_api.cpp 
aointegrals_api.cpp chronusq_api.cpp classtools_api.cpp fileio_api.cpp
tools_api.cpp basisset_api.cpp cerr_api.cpp)

target_link_libraries(pythonapi LINK_PUBLIC basisset classtools controls
mointegrals molecule tools sdresponse realtime mollerplesset singleslater 
aointegrals worker quasinewton grid fileio cerr gauinterface ${Boost_LIBRARIES} 
${PYTHON_LIBRARIES} ${LAPACK_LINKER_FLAGS} ${LAPACK_LIBRARIES} 
${BLAS_LINKER_FLAGS} ${BLAS_LIBRARIES} ${LIBINT_LIB} m ${OpenMP_CXX_FLAGS})

