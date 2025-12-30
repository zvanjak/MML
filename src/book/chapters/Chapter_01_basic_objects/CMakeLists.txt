project(Chapter01App VERSION 1.0)

# Collect all source files for this chapter
file(GLOB CHAPTER_SOURCES "*.cpp")

# Create the executable
add_executable(Chapter01App ${CHAPTER_SOURCES})

# Set output directory
set_target_properties(Chapter01App PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/book_chapters
)
