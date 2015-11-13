add_custom_target(linktests ALL
  COMMAND ${CMAKE_COMMAND} -E create_symlink ${PROJECT_BINARY_DIR}/chronusq.py ${PROJECT_SOURCE_DIR}/tests/chronusq.py
)
