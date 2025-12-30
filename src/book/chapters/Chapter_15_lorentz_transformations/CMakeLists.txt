project(Chapter15App VERSION 1.0)

set(CHAPTER_SOURCES
    main.cpp
    chapter_15_coordinator.cpp
    minkowski_spacetime.cpp
    lorentz_demos.cpp
)

add_executable(Chapter15App ${CHAPTER_SOURCES})
set_target_properties(Chapter15App PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/book_chapters)
