project(Chapter07App VERSION 1.0)

set(CHAPTER_SOURCES
    main.cpp
    chapter_07_coordinator.cpp
    double_pendulum.cpp
    spherical_pendulum.cpp
)
add_executable(Chapter07App ${CHAPTER_SOURCES})
set_target_properties(Chapter07App PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/book_chapters)
