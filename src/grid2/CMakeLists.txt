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
add_library(grid2
SHARED
grid.cpp
lebedev.cpp
#lebedev3.cpp
#lebedev5.cpp
#lebedev7.cpp
#lebedev9.cpp
#lebedev11.cpp
#lebedev13.cpp
#lebedev15.cpp
#lebedev17.cpp
#lebedev19.cpp
lebedev21.cpp
lebedev23.cpp
lebedev25.cpp
lebedev27.cpp
lebedev29.cpp
lebedev31.cpp
lebedev35.cpp
lebedev41.cpp
lebedev47.cpp
lebedev53.cpp
lebedev59.cpp
lebedev65.cpp
#lebedev71.cpp
atomic.cpp)

if(APPLE)
target_link_libraries(grid2 LINK_PUBLIC cerr)
target_link_libraries(grid2 LINK_PUBLIC ${EXTERNAL_CQ_DEPS})
endif()

# Dependencies
add_dependencies(grid2 ${EXTERNAL_CQ_PROJECTS})
