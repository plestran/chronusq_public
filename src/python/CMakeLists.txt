#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I ${PYTHON_INCLUDE_DIRS}")


add_library(pythonapi SHARED singleslater_api.cpp molecule_api.cpp 
aointegrals_api.cpp chronusq_api.cpp classtools_api.cpp fileio_api.cpp
tools_api.cpp basisset_api.cpp cerr_api.cpp)

#target_include_directories(pythonapi PUBLIC ${PYTHON_INCLUDE_DIRS})

target_link_libraries(pythonapi LINK_PUBLIC basisset classtools controls
mointegrals molecule tools sdresponse realtime mollerplesset singleslater 
aointegrals worker quasinewton grid fileio cerr gauinterface ${Boost_LIBRARIES} 
${PYTHON_LIBRARIES} ${LA_LINK} ${LIBINT_LIB} m ${OpenMP_CXX_FLAGS}
${HDF5_CXX_LIBRARIES})

set(APIINST ${PROJECT_BINARY_DIR}/src/python)
set(APIPY   ${APIINST}/chronusq.py)
set(PFOLDER  ${PROJECT_SOURCE_DIR}/src/python/src/parse)
set(JFOLDER  ${PROJECT_SOURCE_DIR}/src/python/src/jobs)

add_custom_target(copyapi ALL
  ${CMAKE_COMMAND} -E copy_if_different ${APIPY} ${PROJECT_BINARY_DIR}
  COMMAND ${CMAKE_COMMAND} -E make_directory ${APIINST}/parse
  COMMAND ${CMAKE_COMMAND} -E make_directory ${APIINST}/jobs
  COMMAND ${CMAKE_COMMAND} -E copy_directory ${PFOLDER} ${APIINST}/parse
  COMMAND ${CMAKE_COMMAND} -E copy_directory ${JFOLDER} ${APIINST}/jobs
)
add_dependencies(copyapi pythonapi)
#file(COPY "${PROJECT_BINARY_DIR}/src/python/chronusq.py"
#  DESTINATION "${PROJECT_BINARY_DIR}"
#  FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
#)
#file(COPY "${PROJECT_SOURCE_DIR}/src/python/src/parse"
#  DESTINATION "${PROJECT_BINARY_DIR}/src/python"
#  FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
#)
#file(COPY "${PROJECT_SOURCE_DIR}/src/python/src/jobs"
#  DESTINATION "${PROJECT_BINARY_DIR}/src/python"
#  FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
#)
