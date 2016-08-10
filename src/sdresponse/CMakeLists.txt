add_library(sdresponse STATIC sdresponse_double.cpp sdresponse_double_io.cpp sdresponse_double_misc.cpp
sdresponse_double_prop.cpp sdresponse_double_qnrelated.cpp sdresponse_complex.cpp 
sdresponse_complex_io.cpp sdresponse_complex_misc.cpp sdresponse_complex_prop.cpp 
sdresponse_complex_qnrelated.cpp)
if(APPLE)
target_link_libraries(sdresponse LINK_PUBLIC singleslater) 
target_link_libraries(sdresponse LINK_PUBLIC quasinewton) 
target_link_libraries(sdresponse LINK_PUBLIC mointegrals) 
target_link_libraries(sdresponse LINK_PUBLIC ${EXTERNAL_CQ_DEPS})
endif()

