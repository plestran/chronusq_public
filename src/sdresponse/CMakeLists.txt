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
add_library(sdresponse STATIC sdresponse_double.cpp sdresponse_double_io.cpp sdresponse_double_misc.cpp
sdresponse_double_prop.cpp sdresponse_double_qnrelated.cpp sdresponse_complex.cpp 
sdresponse_complex_io.cpp sdresponse_complex_misc.cpp sdresponse_complex_prop.cpp 
sdresponse_complex_qnrelated.cpp)
if(APPLE)
target_link_libraries(sdresponse LINK_PUBLIC singleslater) 
target_link_libraries(sdresponse LINK_PUBLIC quasinewton) 
target_link_libraries(sdresponse LINK_PUBLIC mointegrals) 
target_link_libraries(sdresponse LINK_PUBLIC ${EXTERNAL_CQ_DEPS})
endif()

# Dependencies
add_dependencies(sdresponse ${EXTERNAL_CQ_PROJECTS})

