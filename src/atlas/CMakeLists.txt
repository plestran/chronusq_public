add_library(atlas SHARED atlas.cpp)
target_link_libraries(atlas LINK_PUBLIC controls)
target_link_libraries(atlas LINK_PUBLIC classtools)
target_link_libraries(atlas LINK_PUBLIC singleslater)
target_link_libraries(atlas LINK_PUBLIC realtime)
target_link_libraries(atlas LINK_PUBLIC sdresponse)
target_link_libraries(atlas LINK_PUBLIC aointegrals)
target_link_libraries(atlas LINK_PUBLIC mointegrals)
target_link_libraries(atlas LINK_PUBLIC ${EXTERNAL_CQ_DEPS})



if(USE_LIBINT)
if(BUILD_LIBINT)
 add_dependencies(atlas libint)
endif()
endif()
