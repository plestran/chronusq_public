add_library(qn SHARED qn_memory.cpp qn_guess.cpp qn_diag.cpp qn_orth.cpp qn_files.cpp)
target_link_libraries(qn LINK_PUBLIC ${EXTERNAL_CQ_DEPS})
