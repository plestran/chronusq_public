add_library(dft SHARED slater.cpp vwn3.cpp vwn5.cpp b88.cpp lyp.cpp pbe.cpp)

if(APPLE)
target_link_libraries(singleslater LINK_PUBLIC ${EXTERNAL_CQ_DEPS})
endif()

