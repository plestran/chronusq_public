add_library(tools STATIC tools.cpp)
if(APPLE)
target_link_libraries(tools LINK_PUBLIC cerr)
target_link_libraries(tools LINK_PUBLIC ${EXTERNAL_CQ_DEPS})
endif()

