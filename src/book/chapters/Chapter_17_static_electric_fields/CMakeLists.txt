project(Chapter17App VERSION 1.0)

file(GLOB CHAPTER_SOURCES "*.cpp")
add_executable(Chapter17App ${CHAPTER_SOURCES})
set_target_properties(Chapter17App PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/book_chapters)
