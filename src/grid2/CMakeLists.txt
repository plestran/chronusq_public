add_library(grid2 SHARED grid.cpp lebedev.cpp lebedev3.cpp lebedev5.cpp
  lebedev7.cpp lebedev9.cpp lebedev11.cpp lebedev13.cpp lebedev15.cpp
  lebedev17.cpp lebedev19.cpp lebedev21.cpp lebedev23.cpp lebedev25.cpp
  lebedev27.cpp lebedev29.cpp lebedev31.cpp lebedev35.cpp lebedev41.cpp
  lebedev47.cpp lebedev53.cpp lebedev59.cpp lebedev65.cpp lebedev71.cpp)

if(APPLE)
target_link_libraries(grid2 LINK_PUBLIC cerr)
target_link_libraries(grid2 LINK_PUBLIC ${EXTERNAL_CQ_DEPS})
endif()

if(USE_LIBINT)
if(BUILD_LIBINT)
 add_dependencies(grid2 libint)
endif()
endif()
