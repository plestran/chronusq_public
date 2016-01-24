add_library(coupledcluster SHARED coupledcluster.cpp)
if(APPLE)
target_link_libraries(coupledcluster LINK_PUBLIC cerr) 
target_link_libraries(coupledcluster LINK_PUBLIC mointegrals) 
target_link_libraries(coupledcluster LINK_PUBLIC ${EXTERNAL_CQ_DEPS}) 
endif()

if(USE_LIBINT)
if(BUILD_LIBINT)
 add_dependencies(coupledcluster libint)
endif()
endif()

