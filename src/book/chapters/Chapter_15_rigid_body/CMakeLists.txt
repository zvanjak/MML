project(Chapter15App VERSION 1.0)

file(GLOB CHAPTER_SOURCES "*.cpp")
add_executable(Chapter15App ${CHAPTER_SOURCES})
set_target_properties(Chapter15App PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/book_chapters)
