add_library(dft SHARED slater.cpp vwn3.cpp vwn5.cpp b88.cpp lyp.cpp)

if(APPLE)
target_link_libraries(singleslater LINK_PUBLIC ${EXTERNAL_CQ_DEPS})
endif()

if(USE_LIBINT)
if(BUILD_LIBINT)
 add_dependencies(singleslater libint)
endif()
endif()
