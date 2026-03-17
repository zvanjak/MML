project(Chapter12App)

# Collect all source files for this chapter
set(SOURCES
    main.cpp
    chapter_12_coordinator.cpp
    giant_carousel.cpp
    coriolis_demos.cpp
)

# Create executable
add_executable(Chapter12App ${SOURCES})

# Set output directory
set_target_properties(Chapter12App PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/book_chapters"
)
