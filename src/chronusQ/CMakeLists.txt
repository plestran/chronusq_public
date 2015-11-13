add_executable(chronusQ chronusQ.cpp)
#add_dependencies(chronusQ lapack)
#target_link_libraries(chronusQ LINK_PUBLIC atlas classtools basisset controls
#mointegrals molecule tools sdresponse realtime mollerplesset singleslater 
#aointegrals worker quasinewton grid fileio cerr gauinterface 
#${EXTERNAL_CQ_DEPS})

target_link_libraries(chronusQ LINK_PUBLIC atlas)

if(USE_LIBINT)
if(BUILD_LIBINT)
 add_dependencies(chronusQ libint)
endif()
endif()

