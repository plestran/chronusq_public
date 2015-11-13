add_library(quasinewton STATIC quasinewton_real_diag.cpp quasinewton_real_guess.cpp
quasinewton_real_proc.cpp quasinewton_real_restart.cpp
quasinewton_real_misc.cpp quasinewton_complex_diag.cpp quasinewton_complex_guess.cpp
 quasinewton_complex_proc.cpp quasinewton_complex_restart.cpp
quasinewton_complex_misc.cpp)
target_link_libraries(quasinewton LINK_PUBLIC sdresponse)
target_link_libraries(quasinewton LINK_PUBLIC ${EXTERNAL_CQ_DEPS})


if(USE_LIBINT)
if(BUILD_LIBINT)
 add_dependencies(quasinewton libint)
endif()
endif()
