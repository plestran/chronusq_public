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
if(USE_LIBINT)
  set(AOINTS_SOURCE 
     aointegrals.cpp 
     aointegrals_builders.cpp 
#    aointegrals_real_contract_uniquegen.cpp 
#    aointegrals_complex_contract_uniquegen.cpp 
#    aointegrals_real_contract_percontract.cpp 
#    aointegrals_real_contract_direct.cpp 
#    aointegrals_complex_contract_percontract.cpp 
#    aointegrals_complex_contract_direct.cpp aointegrals_real_contract_N4.cpp 
#    aointegrals_real_contract_DF.cpp aointegrals_complex_contract_N4.cpp 
#    aointegrals_complex_contract_DF.cpp aointegrals_dkh.cpp 
    aointegrals_denfit.cpp 
    aointegrals_impl.cpp 
    aointegrals_util.cpp
    aointegrals_onee.cpp
    aointegrals_rel.cpp)
else(USE_LIBINT)
  set(AOINTS_SOURCE aointegrals.cpp aointegrals_builders_nolibint.cpp 
    aointegrals_onee.cpp aointegrals_twoe.cpp) 
endif(USE_LIBINT)

add_library(aointegrals SHARED ${AOINTS_SOURCE})
if(APPLE)
target_link_libraries(aointegrals LINK_PUBLIC basisset)
target_link_libraries(aointegrals LINK_PUBLIC fileio)
target_link_libraries(aointegrals LINK_PUBLIC grid)
target_link_libraries(aointegrals LINK_PUBLIC ${EXTERNAL_CQ_DEPS})
endif()

