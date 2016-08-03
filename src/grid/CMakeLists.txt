add_library(grid SHARED grid.cpp)
if(APPLE)
target_link_libraries(grid LINK_PUBLIC cerr)
target_link_libraries(grid LINK_PUBLIC ${EXTERNAL_CQ_DEPS})
endif()

