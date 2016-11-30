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

# Add an external project to build libint in source
ExternalProject_Add(libint
  PREFIX ${PROJECT_BINARY_DIR}/deps/libint
  URL "${PROJECT_SOURCE_DIR}/deps/src/libint-2.2.0-alpha.tgz"
  CONFIGURE_COMMAND ${PROJECT_BINARY_DIR}/deps/libint/src/libint/configure 
    --prefix=${PROJECT_BINARY_DIR}/deps 
    CXX=${CMAKE_CXX_COMPILER} 
    CXXFLAGS=${CMAKE_CXX_FLAGS} 
    --enable-shared
  BUILD_COMMAND make -j2
  BUILD_IN_SOURCE 1
  INSTALL_COMMAND make install
)

# Manually set the variables that would normally get set by PkgConfig
set(LIBINT2_LIBRARY_DIRS ${PROJECT_BINARY_DIR}/deps/lib)
set(LIBINT2_INCLUDE_DIRS ${PROJECT_BINARY_DIR}/deps/include)
set(LIBINT2_INCLUDE_DIRS ${LIBINT2_INCLUDE_DIRS} ${PROJECT_BINARY_DIR}/deps/include/libint2)

if(APPLE)
  set(LIBINT2_LIBRARIES ${LIBINT2_LIBRARY_DIRS}/libint2.dylib)
else()
  set(LIBINT2_LIBRARIES ${LIBINT2_LIBRARY_DIRS}/libint2.so)
endif()

