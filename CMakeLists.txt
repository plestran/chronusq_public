# Configure Project ChronusQ
cmake_minimum_required(VERSION 2.8)
project(ChronusQ CXX C)
set(ChronusQ_VERSION_MAJOR BETA)
set(ChronusQ_VERSION_MINOR 0)

set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake")

# Check for Dependencies
find_package(BLAS REQUIRED)   # Required for now
find_package(LAPACK REQUIRED) # Required for now
find_package(PythonLibs)
find_package(Boost REQUIRED python)
find_package(Eigen3 REQUIRED)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I ${Boost_INCLUDE_DIRS}")
include_directories("${Boost_INCLUDE_DIRS}")
include_directories("${EIGEN3_INCLUDE_DIR}")

option(USE_OMP "Use SMP Parallelization" OFF) # Default serial processing
option(USE_MPI "Use MPI Parallelization" OFF)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")

if(USE_OMP)
  find_package(OpenMP)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif(USE_OMP)
if(USE_MPI)
  find_package(MPI)
endif(USE_MPI)

option(USE_LIBINT "Use LibInt Integral Driver" OFF) # Default use in-house integrals
option(USE_BOOST_TIMER "Use Boost High Performace Timer" OFF) # Default use built in timer

if(USE_LIBINT)
  find_library(LIBINT_LIB int2 ${PROJECT_SOURCE_DIR}/deps/lib)
  include_directories("${PROJECT_SOURCE_DIR}/deps/include")
  include_directories("${PROJECT_SOURCE_DIR}/deps/include/libint2")
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

