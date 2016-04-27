if(USE_LIBINT)
  set(AOINTS_SOURCE aointegrals.cpp aointegrals_builders.cpp 
    aointegrals_real_contract_uniquegen.cpp 
    aointegrals_complex_contract_uniquegen.cpp 
    aointegrals_real_contract_percontract.cpp 
    aointegrals_real_contract_direct.cpp 
    aointegrals_complex_contract_percontract.cpp 
    aointegrals_complex_contract_direct.cpp aointegrals_real_contract_N4.cpp 
    aointegrals_real_contract_DF.cpp aointegrals_complex_contract_N4.cpp 
    aointegrals_complex_contract_DF.cpp aointegrals_dkh.cpp 
    aointegrals_denfit.cpp aointegrals_util.cpp
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

if(USE_LIBINT)
if(BUILD_LIBINT)
 add_dependencies(aointegrals libint)
endif()
endif()
