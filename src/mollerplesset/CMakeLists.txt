add_library(mollerplesset SHARED mollerplesset.cpp)
if(APPLE)
target_link_libraries(mollerplesset LINK_PUBLIC cerr) 
target_link_libraries(mollerplesset LINK_PUBLIC mointegrals) 
target_link_libraries(mollerplesset LINK_PUBLIC ${EXTERNAL_CQ_DEPS}) 
endif()


