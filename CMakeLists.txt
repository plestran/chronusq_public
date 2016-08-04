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

# Build Needed Parts of Boost
ExternalProject_Add(boost
  PREFIX ${PROJECT_BINARY_DIR}/deps/boost
  URL "${PROJECT_SOURCE_DIR}/deps/src/boost_1_59_0.tar.gz"
  CONFIGURE_COMMAND ${PROJECT_BINARY_DIR}/deps/boost/src/boost/bootstrap.sh
  BUILD_COMMAND ${PROJECT_BINARY_DIR}/deps/boost/src/boost/b2
    --prefix=${PROJECT_BINARY_DIR}/deps 
    --with-python 
    --with-system 
    --ignore-site-config cxxflags=${CMAKE_CXX_FLAGS} 
    install
  BUILD_IN_SOURCE 1
  INSTALL_COMMAND echo "Boost Build Sucess!"
)

# More compile flags that Boost throws a fit about 
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
if(CQ_ENABLE_MPI)
  # MPI Compile flags
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")
else()
  # OpenMP flags
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

if(APPLE)
  set(Boost_LIBRARIES ${PROJECT_BINARY_DIR}/deps/lib/libboost_python.dylib ${PROJECT_BINARY_DIR}/deps/lib/libboost_system.dylib)
else()
  set(Boost_LIBRARIES ${PROJECT_BINARY_DIR}/deps/lib/libboost_python.so ${PROJECT_BINARY_DIR}/deps/lib/libboost_system.so)
endif()

ExternalProject_Add(btas
  PREFIX ${PROJECT_BINARY_DIR}/deps/btas
  GIT_REPOSITORY https://github.com/BTAS/BTAS.git
  UPDATE_COMMAND echo 'Skipping BTAS Update'
  PATCH_COMMAND echo 'Skipping BTAS Patch'
  CONFIGURE_COMMAND echo 'Skipping BTAS Configure'
  BUILD_COMMAND echo 'Skipping BTAS Build'
  INSTALL_COMMAND ${RSYNC_EXECUTABLE} -r 
    ${PROJECT_BINARY_DIR}/deps/btas/src/btas/btas 
    ${PROJECT_BINARY_DIR}/deps/include
)

include_directories("${PROJECT_BINARY_DIR}/deps/include")
link_directories("${PROJECT_BINARY_DIR}/deps/lib")

# Set up Libint environment
#if(USE_LIBINT)
#  if(BUILD_LIBINT)
#    ExternalProject_Add(libint
#      PREFIX ${PROJECT_BINARY_DIR}/deps/libint
#      URL "${PROJECT_SOURCE_DIR}/deps/src/libint-2.2.0-alpha.tgz"
#      CONFIGURE_COMMAND ${PROJECT_BINARY_DIR}/deps/libint/src/libint/configure 
#        --prefix=${PROJECT_BINARY_DIR}/deps 
#        CXX=${CMAKE_CXX_COMPILER} 
#        CXXFLAGS=${CMAKE_CXX_FLAGS} 
#        --enable-shared
#      BUILD_COMMAND make -j2
#      BUILD_IN_SOURCE 1
#      INSTALL_COMMAND make install
#    )
#    if(APPLE)
#      set(LIBINT_LIB ${PROJECT_BINARY_DIR}/deps/lib/libint2.dylib)
#    else()
#      set(LIBINT_LIB ${PROJECT_BINARY_DIR}/deps/lib/libint2.so) 
#    endif()
#    include_directories("${PROJECT_BINARY_DIR}/deps/include")
#    include_directories("${PROJECT_BINARY_DIR}/deps/include/libint2")
#  else()
#    find_library(LIBINT_LIB int2 ${PROJECT_SOURCE_DIR}/deps/lib)
#    include_directories("${PROJECT_SOURCE_DIR}/deps/include")
#    include_directories("${PROJECT_SOURCE_DIR}/deps/include/libint2")
#  endif(BUILD_LIBINT)
#endif(USE_LIBINT)

# To find Libint and LibXC
set( ENV{PKG_CONFIG_PATH} "${PROJECT_SOURCE_DIR}/deps/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")

find_package(PkgConfig)

# Find Libint
pkg_check_modules(LIBINT2 libint2)
link_directories("${LIBINT2_LIBRARY_DIRS}")
include_directories("${LIBINT2_INCLUDE_DIRS}")

if(CQ_ENABLE_LIBXC)
  pkg_check_modules(LIBXC libxc)
  link_directories("${LIBXC_LIBRARY_DIRS}")
endif()


# Optinally build BLAS and LAPACK

if(BUILD_LA)
  ExternalProject_Add(lapack
    PREFIX ${PROJECT_BINARY_DIR}/deps/lapack
    URL "http://www.netlib.org/lapack/lapack-3.5.0.tgz"
    CMAKE_ARGS -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER} 
               -DCMAKE_Fortran_FLAGS='-fPIC'
               -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/deps
  )


  set(LOCAL_BLAS   ${PROJECT_BINARY_DIR}/deps/lib/libblas.a  )
  set(LOCAL_LAPACK ${PROJECT_BINARY_DIR}/deps/lib/liblapack.a)

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
message(STATUS "Will using the following Link Line for Linear Algebra Libs: ${LA_LINK}")

set(EXTERNAL_CQ_DEPS ${Boost_LIBRARIES} ${PYTHON_LIBRARIES} ${LA_LINK} 
  ${LIBINT2_LIBRARIES} ${LIBXC_LIBRARIES} m ${OpenMP_CXX_FLAGS} 
  ${HDF5_CXX_LIBRARIES})

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
