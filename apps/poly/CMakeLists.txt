add_executable(poly2d_main poly2d_main.cc)

target_compile_definitions(poly2d_main PUBLIC -DMAIN)

target_link_libraries(poly2d_main
    acp
    object
    root
    poly
    poly2d
)

add_executable(poly3d_main poly3d_main.cc)

target_compile_definitions(poly3d_main PUBLIC -DMAIN)

target_link_libraries(poly3d_main
    acp
    object
    circle
    encasement2d
    critical_points
    dcel
    predicates2d
    root
    root2
    resultant
    resultant_arrangement
    poly
    poly2d
    poly3d
    ppoly_d
)
