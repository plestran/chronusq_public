add_library(mointegrals SHARED mointegrals_double.cpp mointegrals_complex.cpp)
target_link_libraries(mointegrals LINK_PUBLIC ${EXTERNAL_CQ_DEPS})

if(USE_LIBINT)
if(BUILD_LIBINT)
 add_dependencies(mointegrals libint)
endif()
endif()
