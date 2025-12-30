project(Chapter05App VERSION 1.0)

set(CHAPTER_SOURCES
    main.cpp
    chapter_05_coordinator.cpp
    projectile_demos.cpp
)
add_executable(Chapter05App ${CHAPTER_SOURCES})
set_target_properties(Chapter05App PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/book_chapters)
