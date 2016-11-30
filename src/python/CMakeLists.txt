#
# This file is part of the Chronus Quantum (ChronusQ) software package
# 
# Copyright (C) 2014-2016 Li Research Group (University of Washington)
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
# 
# Contact the Developers:
#   E-Mail: xsli@uw.edu
#
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I ${PYTHON_INCLUDE_DIRS}")


add_library(pythonapi SHARED
singleslater_api.cpp
molecule_api.cpp
aointegrals_api.cpp
chronusq_api.cpp
classtools_api.cpp
fileio_api.cpp
tools_api.cpp
basisset_api.cpp
cerr_api.cpp
#sdresponse_api.cpp
#realtime_api.cpp
#numdiff_api.cpp
)

#target_include_directories(pythonapi PUBLIC ${PYTHON_INCLUDE_DIRS})

target_link_libraries(pythonapi LINK_PUBLIC classtools)
target_link_libraries(pythonapi LINK_PUBLIC quantum)
target_link_libraries(pythonapi LINK_PUBLIC wavefunction)
target_link_libraries(pythonapi LINK_PUBLIC singleslater)
target_link_libraries(pythonapi LINK_PUBLIC quantum)
target_link_libraries(pythonapi LINK_PUBLIC grid2)
target_link_libraries(pythonapi LINK_PUBLIC molecule)
target_link_libraries(pythonapi LINK_PUBLIC tools)
target_link_libraries(pythonapi LINK_PUBLIC cerr)
target_link_libraries(pythonapi LINK_PUBLIC gauinterface)
target_link_libraries(pythonapi LINK_PUBLIC controls)
target_link_libraries(pythonapi LINK_PUBLIC quantum)
target_link_libraries(pythonapi LINK_PUBLIC singleslater)
target_link_libraries(pythonapi LINK_PUBLIC quantum)
target_link_libraries(pythonapi LINK_PUBLIC aointegrals)
target_link_libraries(pythonapi LINK_PUBLIC mointegrals)
target_link_libraries(pythonapi LINK_PUBLIC basisset)
target_link_libraries(pythonapi PUBLIC realtime)
target_link_libraries(pythonapi LINK_PUBLIC qn)
target_link_libraries(pythonapi LINK_PUBLIC dft)
target_link_libraries(pythonapi LINK_PUBLIC fileio)
target_link_libraries(pythonapi LINK_PUBLIC ${EXTERNAL_CQ_DEPS})

# Dependencies
add_dependencies(pythonapi ${EXTERNAL_CQ_PROJECTS})

#target_link_libraries(pythonapi LINK_PUBLIC 
#basisset 
#classtools
#controls
##mointegrals
#molecule
#tools
##sdresponse
##realtime
##mollerplesset
#singleslater
#aointegrals
#worker
##quasinewton
##grid
#fileio
#cerr
#gauinterface
#tools
##response
#qn
##numdiff
#quantum
#wavefunction
#grid2
#dft
#${EXTERNAL_CQ_DEPS})

set(APIINST ${PROJECT_BINARY_DIR}/src/python)
set(APIPY   ${APIINST}/chronusq.py)
set(PFOLDER  ${PROJECT_SOURCE_DIR}/src/python/src/parse)
set(JFOLDER  ${PROJECT_SOURCE_DIR}/src/python/src/jobs)

add_custom_target(copyapi ALL
  ${CMAKE_COMMAND} -E copy_if_different ${APIPY} ${PROJECT_BINARY_DIR}
  COMMAND chmod +x ${PROJECT_BINARY_DIR}/chronusq.py
  COMMAND ${CMAKE_COMMAND} -E make_directory ${APIINST}/parse
  COMMAND ${CMAKE_COMMAND} -E make_directory ${APIINST}/jobs
  COMMAND ${CMAKE_COMMAND} -E copy_directory ${PFOLDER} ${APIINST}/parse
  COMMAND ${CMAKE_COMMAND} -E copy_directory ${JFOLDER} ${APIINST}/jobs
)
add_dependencies(copyapi pythonapi)
if(APPLE)
  add_custom_target(setupMacAPI ALL
    ${CMAKE_COMMAND} -E create_symlink ${PROJECT_BINARY_DIR}/src/python/libpythonapi.dylib ${PROJECT_BINARY_DIR}/src/python/libpythonapi.so
    COMMAND install_name_tool -change libboost_python.dylib ${PROJECT_BINARY_DIR}/deps/lib/libboost_python.dylib ${PROJECT_BINARY_DIR}/src/python/libpythonapi.dylib 
  )
  add_dependencies(setupMacAPI pythonapi)
endif()
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
