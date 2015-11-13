# Configure Project ChronusQ
cmake_minimum_required(VERSION 3.0)
project(ChronusQ CXX C)
set(ChronusQ_VERSION_MAJOR BETA)
set(ChronusQ_VERSION_MINOR 0)

# Append local cmake directory to find CMAKE Modules
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
message(STATUS "Using ${CMAKE_MODULE_PATH}")

# Set default options
option(USE_LIBINT "Use LibInt Integral Driver"  ON) # Default Libint
option(BUILD_LA   "Build LAPACK and BLAS locally" OFF)
option(BUILD_LIBINT "Build Libint via cmake" OFF)
option(NCORE_LIBINT "Number of Cores to build Libint" 4)
option(CQ_ENABLE_ATLAS "Whether or not to search for ATLAS" OFF)
option(BUILD_DOC  "Build Documentation" OFF)

# Check for dependencies
find_package(PythonLibs REQUIRED)
find_package(PythonInterp REQUIRED)
find_package(Boost REQUIRED python)
find_package(Eigen3 REQUIRED)
find_package(HDF5 REQUIRED CXX)
find_package(OpenMP)

find_path(BTAS_INCLUDE_DIR btas/btas.h
  HINTS ${BTAS_ROOT} ENV CPATH)
if(BTAS_INCLUDE_DIR)
  message(STATUS "Found BTAS: " ${BTAS_INCLUDE_DIR})
else()
  message(FATAL "Could not find BTAS" )
endif()

# Define include directories
include_directories("${PROJECT_SOURCE_DIR}/include")
include_directories("${PROJECT_BINARY_DIR}/include")
include_directories("${Boost_INCLUDE_DIRS}")
include_directories("${EIGEN3_INCLUDE_DIR}")
include_directories("${BTAS_INCLUDE_DIR}")
include_directories("${HDF5_INCLUDE_DIRS}")
include_directories("${PYTHON_INCLUDE_DIRS}")

# Append flags for compilation
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

# Set up Libint environment
include(ExternalProject)
if(USE_LIBINT)
  if(BUILD_LIBINT)
    ExternalProject_Add(libint
      PREFIX ${PROJECT_BINARY_DIR}/deps/libint
      URL "${PROJECT_SOURCE_DIR}/deps/src/libint-2.1.0-beta2.tgz"
      CONFIGURE_COMMAND ${PROJECT_BINARY_DIR}/deps/libint/src/libint/configure --prefix=${PROJECT_BINARY_DIR}/deps CXX=${CMAKE_CXX_COMPILER} CXXFLAGS=${CMAKE_CXX_FLAGS}
      BUILD_COMMAND make -j${NCORE_LIBINT}
      BUILD_IN_SOURCE 1
      INSTALL_COMMAND make install
    )  
    set(LIBINT_LIB ${PROJECT_BINARY_DIR}/deps/lib/libint2.a)
    include_directories("${PROJECT_BINARY_DIR}/deps/include")
    include_directories("${PROJECT_BINARY_DIR}/deps/include/libint2")
  else()
    find_library(LIBINT_LIB int2 ${PROJECT_SOURCE_DIR}/deps/lib)
    include_directories("${PROJECT_SOURCE_DIR}/deps/include")
    include_directories("${PROJECT_SOURCE_DIR}/deps/include/libint2")
  endif(BUILD_LIBINT)
endif(USE_LIBINT)

# Optinally build BLAS and LAPACK

if(BUILD_LA)
  ExternalProject_Add(lapack
    PREFIX ${PROJECT_BINARY_DIR}/deps/lapack
    URL "http://www.netlib.org/lapack/lapack-3.5.0.tgz"
    CMAKE_ARGS -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER} 
               -DCMAKE_Fortran_FLAGS='-fPIC'
               -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/deps
  )

# find_library(LOCAL_BLAS blas ${PROJECT_BINARY_DIR}/deps/lib)
# find_library(LOCAL_LAPACK lapack ${PROJECT_BINARY_DIR}/deps/lib)
# find_library(LOCAL_MGLIB tmglib ${PROJECT_BINARY_DIR}/deps/lib)

  set(LOCAL_BLAS   ${PROJECT_BINARY_DIR}/deps/lib/libblas.a  )
  set(LOCAL_LAPACK ${PROJECT_BINARY_DIR}/deps/lib/liblapack.a)
#  set(LOCAL_BLAS ${PROJECT_BINARY_DIR}/deps/lib/libblas.a)

  set(LA_LINK ${LOCAL_LAPACK} ${LOCAL_BLAS} gfortran)
else()
  if(CQ_ENABLE_ATLAS)
    set(BLA_VENDOR ATLAS)
  endif(CQ_ENABLE_ATLAS)

  find_package(BLAS REQUIRED)  
  find_package(LAPACK REQUIRED)
  set(LA_LINK ${LAPACK_LINKER_FLAGS} ${LAPACK_LIBRARIES})
  set(LA_LINK ${LA_LINK} ${BLAS_LINKER_FLAGS} ${BLAS_LIBRARIES})
endif()

# Set configure time variables
set(BASIS_PATH "${PROJECT_SOURCE_DIR}/basis")
set(AUTHOR_PATH "${PROJECT_SOURCE_DIR}/AUTHORS")
set(PYTHONAPI_PATH "${PROJECT_BINARY_DIR}/src/python")
set(USERMANUAL "${CMAKE_BINARY_DIR}/doc/ChronusQUserManual")
set(LOGO "${CMAKE_SOURCE_DIR}/doc/chronus_quantum_logo.png")
set(CHRONUSQ_BIB ${CMAKE_SOURCE_DIR}/doc/ChronusQ)

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
add_subdirectory(doc)
add_subdirectory(tests)
