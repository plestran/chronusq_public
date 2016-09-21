add_library(realtime SHARED realtime_double.cpp realtime_complex.cpp 
realtime_double_envelope.cpp realtime_complex_envelope.cpp 
realtime_double_csv.cpp realtime_complex_csv.cpp realtime_impl.cpp)
if(APPLE)
target_link_libraries(realtime LINK_PUBLIC singleslater)
target_link_libraries(realtime LINK_PUBLIC ${EXTERNAL_CQ_DEPS})
endif()

