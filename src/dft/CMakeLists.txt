add_library(dft SHARED slater.cpp vwn3.cpp)

if(APPLE)
target_link_libraries(singleslater LINK_PUBLIC ${EXTERNAL_CQ_DEPS})
endif()

if(USE_LIBINT)
if(BUILD_LIBINT)
 add_dependencies(singleslater libint)
endif()
endif()
