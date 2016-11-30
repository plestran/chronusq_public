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

# Pull the most recent BTAS from GitHub
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

