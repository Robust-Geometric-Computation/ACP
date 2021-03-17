add_library(critical_points critical_points.cc)
add_library(encasement2d encasement2d.cc)
add_library(predicates2d predicates2d.cc)
add_library(plot plot.cc)
add_library(pngwriter pngwriter.cc)

target_link_libraries(encasement2d
  delaunay
)
