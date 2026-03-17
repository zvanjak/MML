project(Chapter22App)

# Collect all source files for this chapter
file(GLOB SOURCES "*.cpp")

# Create executable
add_executable(Chapter22App ${SOURCES})

# Set output directory
set_target_properties(Chapter22App PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/book_chapters"
)
