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

# Optionally build BLAS and LAPACK
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

# Find LAPACK / BLAS
  if(CQ_ENABLE_ATLAS)
    set(BLA_VENDOR ATLAS)
  endif(CQ_ENABLE_ATLAS)

  find_package(BLAS REQUIRED)  
  find_package(LAPACK REQUIRED)
  set(LA_LINK ${LAPACK_LINKER_FLAGS} ${LAPACK_LIBRARIES})
  set(LA_LINK ${LA_LINK} ${BLAS_LINKER_FLAGS} ${BLAS_LIBRARIES})
endif()
message(STATUS "Will using the following Link Line for Linear Algebra Libs: ${LA_LINK}")

