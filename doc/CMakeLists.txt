find_package(LATEX)

configure_file(
  "${PROJECT_SOURCE_DIR}/doc/ChronusQUserManual.tex.in"
  "${PROJECT_BINARY_DIR}/doc/ChronusQUserManual.tex"
  @ONLY
)

add_custom_target(copydoc ALL
    ${CMAKE_COMMAND} -E copy_if_different ${LOGO} ${CMAKE_BINARY_DIR}/doc
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CHRONUSQ_BIB}.bib 
      ${CMAKE_BINARY_DIR}/doc
)

if(PDFLATEX_COMPILER)
  message(STATUS "PDFLatex Compiler Found! Building Documentaion Target")
  add_custom_target(docs
    ${PDFLATEX_COMPILER} ${USERMANUAL} 
    COMMAND ${BIBTEX_COMMAND} ${USERMANUAL}
    COMMAND ${PDFLATEX_COMPILER} ${USERMANUAL} 
    COMMAND ${PDFLATEX_COMPILER} ${USERMANUAL} 
    DEPENDS copydoc
  )
#  add_dependencies(documentation ${USERMANUAL}.tex)
else()
  message(STATUS "PDFLatex Compiler Not Found! No Building Documentaion Target")
endif()
