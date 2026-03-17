project(Chapter17App VERSION 1.0)

set(CHAPTER_SOURCES
    main.cpp
    chapter_17_coordinator.cpp
    minkowski_spacetime.cpp
    lorentz_demos.cpp
)

add_executable(Chapter17App ${CHAPTER_SOURCES})
set_target_properties(Chapter17App PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/book_chapters)
