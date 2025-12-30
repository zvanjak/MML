project(Chapter09App VERSION 1.0)

set(CHAPTER_SOURCES
    main.cpp
    chapter_09_coordinator.cpp
    five_masses.cpp
    solar_system.cpp
)
add_executable(Chapter09App ${CHAPTER_SOURCES})
set_target_properties(Chapter09App PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/book_chapters)
