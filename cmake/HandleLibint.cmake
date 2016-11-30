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

# Build or Find Libint
if(BUILD_LIBINT)
  message(STATUS "Opting to build Libint2 by CMake")
  include(BuildLibint)
else()
  message(STATUS "Opting to find a pre-compiled Libint")
  include(FindLibint)
endif()

# Append to Include / Link dirs
link_directories("${LIBINT2_LIBRARY_DIRS}")
include_directories("${LIBINT2_INCLUDE_DIRS}")
message(STATUS "Libint2 library will be located at ${LIBINT2_LIBRARY_DIRS}")
message(STATUS "Libint2 headers will be located at ${LIBINT2_INCLUDE_DIRS}")
