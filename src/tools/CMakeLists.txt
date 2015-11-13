add_library(tools SHARED tools.cpp)
target_link_libraries(tools LINK_PUBLIC cerr)
target_link_libraries(tools LINK_PUBLIC ${EXTERNAL_CQ_DEPS})

if(USE_LIBINT)
if(BUILD_LIBINT)
 add_dependencies(tools libint)
endif()
endif()
