add_executable(mattest mattest.cpp mattest_complex.cpp mattest_double.cpp)
target_link_libraries(mattest LINK_PUBLIC matrix tools fileio lapack blas)
