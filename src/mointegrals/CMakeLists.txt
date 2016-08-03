add_library(mointegrals SHARED mointegrals_double.cpp mointegrals_complex.cpp)
if(APPLE)
target_link_libraries(mointegrals LINK_PUBLIC ${EXTERNAL_CQ_DEPS})
endif()

