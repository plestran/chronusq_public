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
add_library(realtime SHARED realtime_double.cpp realtime_complex.cpp 
realtime_double_envelope.cpp realtime_complex_envelope.cpp 
realtime_double_csv.cpp realtime_complex_csv.cpp realtime_impl.cpp)
if(APPLE)
target_link_libraries(realtime LINK_PUBLIC singleslater)
target_link_libraries(realtime LINK_PUBLIC ${EXTERNAL_CQ_DEPS})
endif()

# Dependencies
add_dependencies(realtime ${EXTERNAL_CQ_PROJECTS})
