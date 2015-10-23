add_executable(chronusQ chronusQ.cpp)
target_link_libraries(chronusQ LINK_PUBLIC atlas classtools basisset controls
mointegrals molecule tools sdresponse realtime mollerplesset singleslater 
aointegrals worker quasinewton grid fileio cerr gauinterface coupledcluster 
${LAPACK_LINKER_FLAGS} ${LAPACK_LIBRARIES} ${BLAS_LINKER_FLAGS} 
${BLAS_LIBRARIES} ${LIBINT_LIB} m ${OpenMP_CXX_FLAGS} ${HDF5_CXX_LIBRARIES})

