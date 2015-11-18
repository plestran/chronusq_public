add_library(cerr SHARED cerr.cpp)
target_link_libraries(cerr LINK_PUBLIC ${EXTERNAL_CQ_DEPS})

if(USE_LIBINT)
if(BUILD_LIBINT)
 add_dependencies(cerr libint)
endif()
endif()
