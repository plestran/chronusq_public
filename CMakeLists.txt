#
# This file is part of the Chronus Quantum (ChronusQ) software package
# 
# Copyright (C) 2014-2016 Li Research Group (University of Washington)
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
# 
# Contact the Developers:
#   E-Mail: xsli@uw.edu
#
 

# Configure Project ChronusQ
cmake_minimum_required(VERSION 3.0)
project(ChronusQ CXX C)
set(ChronusQ_VERSION_MAJOR BETA2)
set(ChronusQ_VERSION_MINOR 1)

# Get Git revision hash
if(EXISTS ${PROJECT_SOURCE_DIR}/.git)
  find_package(Git)
  if(GIT_FOUND)
    execute_process(
      COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
      OUTPUT_VARIABLE "CQ_GIT_VERSION"
      ERROR_QUIET
      OUTPUT_STRIP_TRAILING_WHITESPACE
    )
  else(GIT_FOUND)
    set(CQ_GIT_VERSION 0)
  endif(GIT_FOUND)
endif()

# Append local cmake directory to find CMAKE Modules
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
message(STATUS "Using ${CMAKE_MODULE_PATH}")

# Set default options
option(USE_LIBINT "Use LibInt Integral Driver"  ON) # Default Libint
option(BUILD_LA   "Build LAPACK and BLAS locally" OFF)
option(BUILD_LIBINT "Build Libint via cmake" ON)
option(CQ_ENABLE_ATLAS "Whether or not to search for ATLAS" OFF)
option(CQ_ENABLE_MPI "Enable MPI Parallelism" OFF)
option(CQ_ENABLE_LIBXC "Enable LibXC for DFT calculations" OFF)

# Check for C++11
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
  message(FATAL "${CMAKE_CXX_COMPILER} does not seem to support C++11!")
endif()

# Check for dependencies
find_package(PythonLibs 2.7 REQUIRED)
find_package(PythonInterp 2.7 REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(HDF5 REQUIRED CXX)
find_program(RSYNC_EXECUTABLE rsync)

if(CQ_ENABLE_MPI)
  find_package(MPI)
  if(NOT MPI_CXX_FOUND)
    message(FATAL "Suitable CXX binding for MPI not found")
  endif()
else()
  find_package(OpenMP)
endif()


# Define include directories
include_directories("${PROJECT_SOURCE_DIR}/include")
include_directories("${PROJECT_BINARY_DIR}/include")
include_directories("${EIGEN3_INCLUDE_DIR}")
include_directories("${HDF5_INCLUDE_DIRS}")
include_directories("${PYTHON_INCLUDE_DIRS}")
if(CQ_ENABLE_MPI)
  include_directories("${MPI_CXX_INCLUDE_PATH}")
endif()


if(NOT CMAKE_CXX_COMPILER_ID MATCHES "PGI")
  # Append flags for compilation
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-declarations")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-narrowing ")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-strict-aliasing")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-variable")
  
  if(NOT CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-local-typedefs")
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-return-type")
  endif()
endif()

# Vectorization
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftree-vectorize")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse2")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffast-math")

include(ExternalProject)

# Boost and BTAS Always built
include(BuildBoost)
include(BuildBTAS)


# More compile flags that Boost throws a fit about 
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
if(CQ_ENABLE_MPI)
  # MPI Compile flags
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")
else()
  # OpenMP flags
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

include_directories("${PROJECT_BINARY_DIR}/deps/include")
link_directories("${PROJECT_BINARY_DIR}/deps/lib")


include(FindBuiltLibs)
include(FindLA)

set(EXTERNAL_CQ_DEPS ${Boost_LIBRARIES} ${PYTHON_LIBRARIES} ${LA_LINK} 
  ${LIBINT2_LIBRARIES} ${LIBXC_LIBRARIES} m ${OpenMP_CXX_FLAGS} 
  ${HDF5_CXX_LIBRARIES})

set(EXTERNAL_CQ_PROJECTS boost btas)
if(BUILD_LIBINT)
  set(EXTERNAL_CQ_PROJECTS ${EXTERNAL_CQ_PROJECTS} libint)
endif()

if(CQ_ENABLE_MPI)
  set(EXTERNAL_CQ_DEPS ${EXTERNAL_CQ_DEPS} ${MPI_CXX_LIBRARIES})
endif()

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
add_subdirectory(tests)


message(STATUS "CXXFLAGS: ${CMAKE_CXX_FLAGS}")
