add_library(worker SHARED worker.cpp)
if(APPLE)
target_link_libraries(worker LINK_PUBLIC ${EXTERNAL_CQ_DEPS})
endif()
