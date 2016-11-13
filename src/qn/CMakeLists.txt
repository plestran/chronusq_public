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
add_library(qn SHARED qn_memory.cpp qn_guess.cpp qn_diag.cpp qn_orth.cpp qn_files.cpp qn_special.cpp)
if(APPLE)
target_link_libraries(qn LINK_PUBLIC classtools)
target_link_libraries(qn LINK_PUBLIC tools)
target_link_libraries(qn LINK_PUBLIC cerr)
endif()

target_link_libraries(qn LINK_PUBLIC ${EXTERNAL_CQ_DEPS})
