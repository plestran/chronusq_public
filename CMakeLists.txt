cmake_minimum_required(VERSION 3.0)
project(ChronusQ CXX Fortran)
set(ChronusQ_VERSION_MAJOR 1)
set(ChronusQ_VERSION_MINOR 0)

configure_file(
  "${PROJECT_SOURCE_DIR}/include/config.h.in"
  "${PROJECT_BINARY_DIR}/include/config.h"
)

include_directories("${PROJECT_SOURCE_DIR}/include")
include_directories("${PROJECT_BINARY_DIR}/include")

add_subdirectory(src)

