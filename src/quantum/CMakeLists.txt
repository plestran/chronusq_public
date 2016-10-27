add_library(quantum STATIC 
#quantum_copy.cpp 
quantum_impl.cpp)

if(APPLE)
target_link_libraries(quantum LINK_PUBLIC singleslater)
target_link_libraries(quantum LINK_PUBLIC ${EXTERNAL_CQ_DEPS})
endif()

