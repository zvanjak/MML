project(Chapter11App)

# Collect all source files for this chapter
set(SOURCES
    main.cpp
    chapter_11_coordinator.cpp
    giant_carousel.cpp
    coriolis_demos.cpp
)

# Create executable
add_executable(Chapter11App ${SOURCES})

# Set output directory
set_target_properties(Chapter11App PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/book_chapters"
)
