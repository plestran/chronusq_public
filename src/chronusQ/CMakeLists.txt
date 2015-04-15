add_executable(chronusQ chronusQ.cpp)
target_link_libraries(chronusQ LINK_PUBLIC atlas basisset classtools controls
mointegrals molecule tools sdresponse singleslater aointegrals worker matrix
fileio lapack blas ${LIBINT_LIB})

