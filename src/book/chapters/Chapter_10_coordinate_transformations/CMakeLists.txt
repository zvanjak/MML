project(Chapter10App VERSION 1.0)

set(CHAPTER_SOURCES
    main.cpp
    chapter_10_coordinator.cpp
    oblique_transforms.cpp
    spherical_transforms.cpp
    cylindrical_transforms.cpp
    basis_vectors.cpp
)

add_executable(Chapter10App ${CHAPTER_SOURCES})
set_target_properties(Chapter10App PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/book_chapters)
