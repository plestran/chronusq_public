add_library(basisset SHARED basisset_gen.cpp basisset_map.cpp basisset_misc.cpp
basisset_print.cpp basisset_read.cpp)
target_link_libraries(basisset LINK_PUBLIC tools)
target_link_libraries(basisset LINK_PUBLIC cerr)
target_link_libraries(basisset LINK_PUBLIC ${EXTERNAL_CQ_DEPS})


if(USE_LIBINT)
if(BUILD_LIBINT)
 add_dependencies(basisset libint)
endif()
endif()
