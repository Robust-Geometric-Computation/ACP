add_library(poly poly.cc)
add_library(ppoly ppoly.cc)
add_library(ppoly_d ppoly_d.cc)
add_library(root root.cc)
add_library(root2 root2.cc)

target_link_libraries(poly root)
target_link_libraries(root poly)

add_library(poly2d poly2d.cc)
add_library(poly3d poly3d.cc)

target_link_libraries(poly2d
  poly
  root
)

target_link_libraries(poly3d
  poly
  poly2d
  root
)
