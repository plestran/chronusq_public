add_executable(chronusQ chronusQ.cpp)
target_link_libraries(chronusQ LINK_PUBLIC atlas basisset classtools controls
mointegrals molecule tools sdresponse singleslater aointegrals worker
fileio lapack blas ${LIBINT_LIB} m gmpxx gmp ${OpenMP_CXX_FLAGS})

