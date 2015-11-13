add_library(controls SHARED controls.cpp)
target_link_libraries(controls LINK_PUBLIC tools)
target_link_libraries(controls LINK_PUBLIC ${EXTERNAL_CQ_DEPS})

if(USE_LIBINT)
if(BUILD_LIBINT)
 add_dependencies(controls libint)
endif()
endif()
