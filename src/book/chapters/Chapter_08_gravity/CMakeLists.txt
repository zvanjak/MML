project(Chapter08App VERSION 1.0)

set(CHAPTER_SOURCES
    main.cpp
    chapter_08_coordinator.cpp
    two_body_2d.cpp
    two_body_3d.cpp
    field_operations.cpp
)

add_executable(Chapter08App ${CHAPTER_SOURCES})
set_target_properties(Chapter08App PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/book_chapters)
