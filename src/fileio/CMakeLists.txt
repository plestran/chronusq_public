add_library(fileio SHARED fileio.cpp)
if(APPLE)
target_link_libraries(fileio LINK_PUBLIC cerr)
target_link_libraries(fileio LINK_PUBLIC ${EXTERNAL_CQ_DEPS})
endif()

