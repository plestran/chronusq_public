add_executable(chronusQ chronusQ.cpp)
target_link_libraries(chronusQ LINK_PUBLIC atlas classtools basisset controls
mointegrals molecule tools sdresponse mollerplesset singleslater aointegrals worker
quasinewton fileio cerr gauinterface lapack blas ${LIBINT_LIB} m ${OpenMP_CXX_FLAGS})

