add_executable(descartes descartes.cc)

target_link_libraries(descartes
    acp
    object
    root
    poly
)
