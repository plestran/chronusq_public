add_library(controls SHARED controls.cpp)
if(APPLE)
target_link_libraries(controls LINK_PUBLIC tools)
target_link_libraries(controls LINK_PUBLIC ${EXTERNAL_CQ_DEPS})
endif()

