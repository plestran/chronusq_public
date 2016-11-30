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

# Build Needed Parts of Boost
# FIXME: Should move this to the modular boost from the GitHub Repo
ExternalProject_Add(boost
  PREFIX ${PROJECT_BINARY_DIR}/deps/boost
  URL "${PROJECT_SOURCE_DIR}/deps/src/boost_1_59_0.tar.gz"
  CONFIGURE_COMMAND ${PROJECT_BINARY_DIR}/deps/boost/src/boost/bootstrap.sh
  BUILD_COMMAND ${PROJECT_BINARY_DIR}/deps/boost/src/boost/b2
    --prefix=${PROJECT_BINARY_DIR}/deps 
    --with-python 
    --with-system 
    --with-math
    --ignore-site-config cxxflags=${CMAKE_CXX_FLAGS} 
    install
  BUILD_IN_SOURCE 1
  INSTALL_COMMAND echo "Boost Build Sucess!"
)

# Set vars to define linked boost libs
if(APPLE)
  set(Boost_LIBRARIES ${PROJECT_BINARY_DIR}/deps/lib/libboost_python.dylib ${PROJECT_BINARY_DIR}/deps/lib/libboost_system.dylib)
else()
  set(Boost_LIBRARIES 
      ${PROJECT_BINARY_DIR}/deps/lib/libboost_python.so 
      ${PROJECT_BINARY_DIR}/deps/lib/libboost_system.so
      ${PROJECT_BINARY_DIR}/deps/lib/libboost_math_tr1.so
      ${PROJECT_BINARY_DIR}/deps/lib/libboost_math_tr1f.so
      ${PROJECT_BINARY_DIR}/deps/lib/libboost_math_tr1l.so
      ${PROJECT_BINARY_DIR}/deps/lib/libboost_math_c99.so
      ${PROJECT_BINARY_DIR}/deps/lib/libboost_math_c99f.so
      ${PROJECT_BINARY_DIR}/deps/lib/libboost_math_c99l.so)
endif()
