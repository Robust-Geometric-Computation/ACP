add_library(resultant resultant.cc)
add_library(resultant_arrangement resultant_arrangement.cc)

target_link_libraries(resultant
  predicates2d
)

target_link_libraries(resultant_arrangement
  resultant
)
