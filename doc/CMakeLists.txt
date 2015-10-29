find_package(LATEX)

set(USERMANUAL "${CMAKE_BINARY_DIR}/doc/ChronusQUserManual")
set(LOGO "${CMAKE_SOURCE_DIR}/doc/chronus_quantum_logo.png")

if(PDFLATEX_COMPILER)
  message(STATUS "PDFLatex Compiler Found! Building Documentaion")
  configure_file(
    "${PROJECT_SOURCE_DIR}/doc/ChronusQUserManual.tex.in"
    "${PROJECT_BINARY_DIR}/doc/ChronusQUserManual.tex"
  )
  add_custom_target(documentation ALL
    ${CMAKE_COMMAND} -E copy_if_different ${LOGO} ${CMAKE_BINARY_DIR}/doc
    COMMAND ${PDFLATEX_COMPILER} ${USERMANUAL} 
    COMMAND ${PDFLATEX_COMPILER} ${USERMANUAL} 
  )
#  add_dependencies(documentation ${USERMANUAL}.tex)
else()
  message(STATUS "PDFLatex Compiler Not Found! Not Building Documentaion")
endif()
