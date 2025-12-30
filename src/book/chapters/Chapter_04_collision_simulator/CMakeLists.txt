project(Chapter04App VERSION 1.0)

file(GLOB CHAPTER_SOURCES "*.cpp")
add_executable(Chapter04App ${CHAPTER_SOURCES})
set_target_properties(Chapter04App PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/book_chapters)
