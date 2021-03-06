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
add_executable(gentests gentests2.cpp)
target_link_libraries(gentests LINK_PUBLIC classtools)
target_link_libraries(gentests LINK_PUBLIC quantum)
target_link_libraries(gentests LINK_PUBLIC wavefunction)
target_link_libraries(gentests LINK_PUBLIC grid2)
target_link_libraries(gentests LINK_PUBLIC singleslater)
target_link_libraries(gentests LINK_PUBLIC quantum)
#target_link_libraries(gentests LINK_PUBLIC grid)
target_link_libraries(gentests LINK_PUBLIC grid2)
target_link_libraries(gentests LINK_PUBLIC molecule)
target_link_libraries(gentests LINK_PUBLIC tools)
target_link_libraries(gentests LINK_PUBLIC cerr)
target_link_libraries(gentests LINK_PUBLIC gauinterface)
target_link_libraries(gentests LINK_PUBLIC controls)
target_link_libraries(gentests LINK_PUBLIC quantum)
target_link_libraries(gentests LINK_PUBLIC singleslater)
target_link_libraries(gentests LINK_PUBLIC grid2)
target_link_libraries(gentests LINK_PUBLIC wavefunction)
target_link_libraries(gentests LINK_PUBLIC quantum)
target_link_libraries(gentests LINK_PUBLIC aointegrals)
target_link_libraries(gentests LINK_PUBLIC mointegrals)
target_link_libraries(gentests LINK_PUBLIC basisset)
#target_link_libraries(gentests LINK_PUBLIC response)
target_link_libraries(gentests PUBLIC realtime)
target_link_libraries(gentests LINK_PUBLIC qn)
target_link_libraries(gentests LINK_PUBLIC dft)
target_link_libraries(gentests LINK_PUBLIC fileio)
target_link_libraries(gentests LINK_PUBLIC ${EXTERNAL_CQ_DEPS})

# Dependencies
add_dependencies(gentests ${EXTERNAL_CQ_PROJECTS})
