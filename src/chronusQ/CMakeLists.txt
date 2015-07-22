add_executable(chronusQ chronusQ.cpp)
target_link_libraries(chronusQ LINK_PUBLIC atlas classtools basisset controls
mointegrals molecule tools sdresponse singleslater aointegrals worker
quasinewton grid fileio cerr gauinterface lapack blas ${LIBINT_LIB} m ${OpenMP_CXX_FLAGS})

