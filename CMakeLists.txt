# Configure Project ChronusQ
cmake_minimum_required(VERSION 2.8)
project(ChronusQ C CXX Fortran)
set(ChronusQ_VERSION_MAJOR 1)
set(ChronusQ_VERSION_MINOR 0)

# Check for Dependencies
find_package(BLAS REQUIRED)   # Required for now
find_package(LAPACK REQUIRED) # Required for now

option(USE_OMP "Use SMP Parallelization" OFF) # Default serial processing
option(USE_MPI "Use MPI Parallelization" OFF)

if(USE_OMP)
  find_package(OpenMP)
endif(USE_OMP)
if(USE_MPI)
  find_package(MPI)
endif(USE_MPI)

option(USE_LIBINT "Use LibInt Integral Driver" OFF) # Default use in-house integrals
option(USE_BOOST_TIMER "Use Boost High Performace Timer" OFF) # Default use built in timer

if(USE_LIBINT)
  find_library(LIBINT_LIB int2 /usr/local/libint/2.1.0-beta/lib)
  include_directories("/usr/local/libint/2.1.0-beta/include")
  include_directories("/usr/local/libint/2.1.0-beta/include/libint2")
  include_directories("/usr/include/eigen3")
endif(USE_LIBINT)


set(BASIS_PATH "${PROJECT_SOURCE_DIR}/basis")

# Set up config.h with configure time settings
configure_file(
  "${PROJECT_SOURCE_DIR}/include/config.h.in"
  "${PROJECT_BINARY_DIR}/include/config_chronusq.h"
)

# Define include Dir
include_directories("${PROJECT_SOURCE_DIR}/include")
include_directories("${PROJECT_BINARY_DIR}/include")

# Recurse into src
add_subdirectory(src)

