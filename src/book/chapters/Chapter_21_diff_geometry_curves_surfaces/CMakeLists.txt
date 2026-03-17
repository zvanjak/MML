project(Chapter21App)

# Collect all source files for this chapter
set(SOURCES
    main.cpp
    chapter_21_coordinator.cpp
    curves_demo.cpp
    surfaces_demo.cpp
)

# Create executable
add_executable(Chapter21App ${SOURCES})

# Set output directory
set_target_properties(Chapter21App PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/book_chapters"
)
