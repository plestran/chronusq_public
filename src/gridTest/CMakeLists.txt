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
add_executable(GRIDTest test.cxx)
target_link_libraries(GRIDTest PUBLIC classtools)
target_link_libraries(GRIDTest PUBLIC controls)
target_link_libraries(GRIDTest PUBLIC singleslater)
target_link_libraries(GRIDTest PUBLIC aointegrals)
target_link_libraries(GRIDTest PUBLIC mointegrals)
target_link_libraries(GRIDTest PUBLIC basisset)
target_link_libraries(GRIDTest PUBLIC response)
target_link_libraries(GRIDTest PUBLIC fileio)
target_link_libraries(GRIDTest PUBLIC ${EXTERNAL_CQ_DEPS})

# Dependencies
add_dependencies(GRIDTest ${EXTERNAL_CQ_PROJECTS})
