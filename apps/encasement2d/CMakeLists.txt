add_executable(encasement_plot encasement_plot.cc)

target_link_libraries(encasement_plot
    acp
    object
    delaunay
    encasement2d
    encasement_utils
    circle
    critical_points
    dcel
    predicates2d
    plot
    pngwriter
    root
    root2
    resultant_arrangement
    poly
    poly2d
    ppoly_d
)

target_link_libraries(encasement_plot
    ${OPENGL_LIBRARIES}
    ${GLUT_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${FREETYPE_LIBRARIES}
    ${PNG_LIBRARIES}
    ${ZLIB_LIBRARIES}
)

add_executable(encasement_test encasement_test.cc)

target_link_libraries(encasement_test
    acp
    object
    delaunay
    encasement2d
    circle
    critical_points
    dcel
    predicates2d
    root
    root2
    resultant_arrangement
    poly
    poly2d
    ppoly_d
)

target_link_libraries(encasement_test
    ${LAPACK_LIBRARIES}
)
