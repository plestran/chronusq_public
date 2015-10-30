find_package(LATEX)


if(PDFLATEX_COMPILER)
  message(STATUS "PDFLatex Compiler Found! Building Documentaion")
  add_custom_target(documentation ALL
    ${CMAKE_COMMAND} -E copy_if_different ${LOGO} ${CMAKE_BINARY_DIR}/doc
    COMMAND ${PDFLATEX_COMPILER} ${USERMANUAL} 
    COMMAND ${PDFLATEX_COMPILER} ${USERMANUAL} 
  )
#  add_dependencies(documentation ${USERMANUAL}.tex)
else()
  message(STATUS "PDFLatex Compiler Not Found! Not Building Documentaion")
endif()
