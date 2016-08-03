add_library(basisset SHARED basisset_gen.cpp basisset_map.cpp basisset_misc.cpp
basisset_print.cpp basisset_read.cpp)
if(APPLE)
target_link_libraries(basisset LINK_PUBLIC tools)
target_link_libraries(basisset LINK_PUBLIC cerr)
target_link_libraries(basisset LINK_PUBLIC ${EXTERNAL_CQ_DEPS})
endif()


