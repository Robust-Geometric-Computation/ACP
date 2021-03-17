## "surface" is the name of the executable

add_executable(display_mesh_folder_vtk display_mesh_folder_vtk.cc)
target_link_libraries(display_mesh_folder_vtk ${VTK_LIBRARIES})

add_executable(surface_visualization_vtk surface_visualization_vtk.cc)
target_link_libraries(surface_visualization_vtk
    acp
    object
    encasement_utils
    poly3d
    ${VTK_LIBRARIES}
    ${GFLAGS_LIBRARIES}
)
add_custom_command(TARGET surface_visualization_vtk POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                       ${CMAKE_SOURCE_DIR}/extras/gradient.jpg
                       $<TARGET_FILE_DIR:surface_visualization_vtk>
                   )
add_custom_command(TARGET surface_visualization_vtk POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                       ${CMAKE_SOURCE_DIR}/extras/paraview.jpg
                       $<TARGET_FILE_DIR:surface_visualization_vtk>
                   )
set(SOURCE_FILES
    viz3d.cc)

add_executable(surface ${SOURCE_FILES} ${GLAD})

# acp targets
target_link_libraries(surface
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
    visualization
)

add_executable(meshviz3d meshviz3d.cc ${GLAD})

# acp targets
target_link_libraries(meshviz3d
    acp
    object
    delaunay
    encasement2d
    encasement3d
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
    poly3d
    ppoly_d
    visualization
    )

  target_link_libraries(meshviz3d
        ${LAPACK_LIBRARIES}
        )



# platform specific targets
if(APPLE)
  target_link_libraries (surface
          glfw
          ${OPENGL_LIBRARIES}
          ${GLEW_LIBRARY}
          ${GMP_LIBRARIES}
          ${COCOA_LIBRARY} ${COREVID_LIBRARY} ${IOKIT_LIBRARY})
  target_link_libraries (meshviz3d
          glfw
          ${OPENGL_LIBRARIES}
          ${GLEW_LIBRARY}
          ${GMP_LIBRARIES}
          ${COCOA_LIBRARY} ${COREVID_LIBRARY} ${IOKIT_LIBRARY})
endif(APPLE)

if (LINUX)
  target_link_libraries (surface
          ${OPENGL_LIBRARIES}
          glfw
          ${GLEW_LIBRARIES}
          Xxf86vm pthread Xrandr Xinerama Xi Xcursor)
  target_link_libraries (meshviz3d
          ${OPENGL_LIBRARIES}
          glfw
          ${GLEW_LIBRARIES}
          Xxf86vm pthread Xrandr Xinerama Xi Xcursor)
endif (LINUX)

