add_library(fileio SHARED fileio.cpp)
target_link_libraries(fileio LINK_PUBLIC cerr)
target_link_libraries(fileio LINK_PUBLIC ${EXTERNAL_CQ_DEPS})

if(USE_LIBINT)
if(BUILD_LIBINT)
 add_dependencies(fileio libint)
endif()
endif()
