set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I /usr/include/python2.6")
add_library(pythonapi SHARED singleslater_api.cpp)
target_link_libraries(pythonapi LINK_PUBLIC atlas classtools basisset controls
mointegrals molecule tools sdresponse realtime mollerplesset singleslater aointegrals worker
quasinewton grid fileio cerr gauinterface m boost_python python2.6)

