add_library(acp acp.cc)
target_link_libraries(acp ${MPFR_LIBRARIES} ${GMP_LIBRARIES} Threads::Threads)

add_library(object object.cc)
