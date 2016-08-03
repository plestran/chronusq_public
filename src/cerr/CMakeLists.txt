add_library(cerr STATIC cerr.cpp)
if(APPLE)
target_link_libraries(cerr LINK_PUBLIC ${EXTERNAL_CQ_DEPS})
target_link_libraries(cerr LINK_PUBLIC tools)
endif()

