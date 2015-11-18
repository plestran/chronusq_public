add_library(grid SHARED grid.cpp)
target_link_libraries(grid LINK_PUBLIC cerr)
target_link_libraries(grid LINK_PUBLIC ${EXTERNAL_CQ_DEPS})

if(USE_LIBINT)
if(BUILD_LIBINT)
 add_dependencies(grid libint)
endif()
endif()
