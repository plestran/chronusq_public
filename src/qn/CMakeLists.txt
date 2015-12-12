add_library(qn SHARED qn_memory.cpp qn_guess.cpp)
target_link_libraries(qn LINK_PUBLIC ${EXTERNAL_CQ_DEPS})
