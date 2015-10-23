#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I ${PYTHON_INCLUDE_DIRS}")


add_library(pythonapi SHARED singleslater_api.cpp molecule_api.cpp 
aointegrals_api.cpp chronusq_api.cpp classtools_api.cpp fileio_api.cpp
tools_api.cpp basisset_api.cpp cerr_api.cpp)

target_include_directories(pythonapi PUBLIC ${PYTHON_INCLUDE_DIRS})

target_link_libraries(pythonapi LINK_PUBLIC basisset classtools controls
mointegrals molecule tools sdresponse realtime mollerplesset singleslater 
aointegrals worker quasinewton grid fileio cerr gauinterface ${Boost_LIBRARIES} 
${PYTHON_LIBRARIES} ${LAPACK_LINKER_FLAGS} ${LAPACK_LIBRARIES} 
${BLAS_LINKER_FLAGS} ${BLAS_LIBRARIES} ${LIBINT_LIB} m ${OpenMP_CXX_FLAGS}
${HDF5_CXX_LIBRARIES})

file(COPY "${PROJECT_BINARY_DIR}/src/python/chronusq.py"
  DESTINATION "${PROJECT_BINARY_DIR}"
  FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
)
file(COPY "${PROJECT_SOURCE_DIR}/src/python/src/parse"
  DESTINATION "${PROJECT_BINARY_DIR}/src/python"
  FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
)
file(COPY "${PROJECT_SOURCE_DIR}/src/python/src/jobs"
  DESTINATION "${PROJECT_BINARY_DIR}/src/python"
  FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
)
