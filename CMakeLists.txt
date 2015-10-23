# Configure Project ChronusQ
cmake_minimum_required(VERSION 2.8)
project(ChronusQ CXX C)
set(ChronusQ_VERSION_MAJOR BETA)
set(ChronusQ_VERSION_MINOR 0)

# Append local cmake directory to find CMAKE Modules
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake")

# Set default options
option(USE_LIBINT "Use LibInt Integral Driver"  ON) # Default Libint

# Check for dependencies
find_package(BLAS REQUIRED)  
find_package(LAPACK REQUIRED)
find_package(PythonLibs REQUIRED)
find_package(PythonInterp REQUIRED)
find_package(Boost REQUIRED python)
find_package(Eigen3 REQUIRED)
find_package(HDF5 REQUIRED)
find_package(OpenMP)

# Define include directories
include_directories("${PROJECT_SOURCE_DIR}/include")
include_directories("${PROJECT_BINARY_DIR}/include")
include_directories("${Boost_INCLUDE_DIRS}")
include_directories("${EIGEN3_INCLUDE_DIR}")
include_directories("${HDF5_INCLUDE_DIRS}")

# Append flags for compilation
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

# Set up Libint environment
if(USE_LIBINT)
  find_library(LIBINT_LIB int2 ${PROJECT_SOURCE_DIR}/deps/lib)
  include_directories("${PROJECT_SOURCE_DIR}/deps/include")
  include_directories("${PROJECT_SOURCE_DIR}/deps/include/libint2")
endif(USE_LIBINT)

# Set configure time variables
set(BASIS_PATH "${PROJECT_SOURCE_DIR}/basis")
set(AUTHOR_PATH "${PROJECT_SOURCE_DIR}/AUTHORS")
set(PYTHONAPI_PATH "${PROJECT_BINARY_DIR}/src/python")

# Set up config.h with configure time settings
configure_file(
  "${PROJECT_SOURCE_DIR}/include/config.h.in"
  "${PROJECT_BINARY_DIR}/include/config_chronusq.h"
)
configure_file(
  "${PROJECT_SOURCE_DIR}/src/python/src/chronusq.py.in"
  "${PROJECT_BINARY_DIR}/src/python/chronusq.py"
)


# Recurse into src
add_subdirectory(src)
