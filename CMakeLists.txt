cmake_minimum_required(VERSION 2.8)
project(ChronusQ C CXX Fortran)
set(ChronusQ_VERSION_MAJOR 1)
set(ChronusQ_VERSION_MINOR 0)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(OpenMP)
find_package(MPI)
configure_file(
  "${PROJECT_SOURCE_DIR}/include/config.h.in"
  "${PROJECT_BINARY_DIR}/include/config.h"
)

include_directories("${PROJECT_SOURCE_DIR}/include")
include_directories("${PROJECT_BINARY_DIR}/include")

add_subdirectory(src)

