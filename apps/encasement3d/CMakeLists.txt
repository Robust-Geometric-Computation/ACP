add_executable(encasement3d_getRoots_test encasement3d_getRoots_test.cc)

target_link_libraries(encasement3d_getRoots_test
    acp
    object
    encasement_utils
    poly3d
)

target_link_libraries(encasement3d_getRoots_test
    ${LAPACK_LIBRARIES}
)

add_executable(encasement3d_generate encasement3d_generate.cc)

target_link_libraries(encasement3d_generate
        acp
        object
        encasement_utils
        poly3d
        )

target_link_libraries(encasement3d_generate
        ${LAPACK_LIBRARIES}
        )

add_executable(encasement3d_test encasement3d_test.cc)

target_link_libraries(encasement3d_test
        acp
        object
        encasement_utils
        poly3d
        )

target_link_libraries(encasement3d_test
        ${LAPACK_LIBRARIES}
        )

add_executable(splitting_plane_test splitting_plane_test.cc)

target_link_libraries(splitting_plane_test
        acp
        object
        delaunay
        encasement2d
        encasement_utils
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
        poly3d
        )

target_link_libraries(splitting_plane_test
        ${LAPACK_LIBRARIES}
        )
