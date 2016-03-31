add_library(singleslater SHARED singleslater_real_fock.cpp singleslater_real_guess.cpp 
singleslater_real_mpi.cpp singleslater_real_scf.cpp 
singleslater_real_misc.cpp singleslater_real_diis.cpp singleslater_complex_misc.cpp 
singleslater_complex_guess.cpp singleslater_complex_diis.cpp singleslater_complex_scf.cpp
singleslater_real_print.cpp singleslater_complex_print.cpp
singleslater_real_dft.cpp singleslater_complex_dft.cpp)

if(APPLE)
target_link_libraries(singleslater LINK_PUBLIC basisset)
target_link_libraries(singleslater LINK_PUBLIC controls)
target_link_libraries(singleslater LINK_PUBLIC molecule)
target_link_libraries(singleslater LINK_PUBLIC grid)
target_link_libraries(singleslater LINK_PUBLIC quantum)
target_link_libraries(singleslater LINK_PUBLIC aointegrals)
target_link_libraries(singleslater LINK_PUBLIC gauinterface)
target_link_libraries(singleslater LINK_PUBLIC ${EXTERNAL_CQ_DEPS})
endif()

if(USE_LIBINT)
if(BUILD_LIBINT)
 add_dependencies(singleslater libint)
endif()
endif()
