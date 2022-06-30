cmake_minimum_required(VERSION 3.0.0)

set (CMAKE_CXX_STANDARD 20)

project(MML_Demo VERSION 0.1.0)

set(MML_DEMO_SOURCES src/main.cpp 
                     src/demo_vector.cpp )
add_executable(MML_Demo ${MML_DEMO_SOURCES})

add_compile_definitions(MML_USE_SINGLE_HEADER)
target_include_directories(MML_Demo PRIVATE include)

set(CATCH_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/catch)
add_library(Catch INTERFACE)
target_include_directories(Catch INTERFACE ${CATCH_INCLUDE_DIR})

file(GLOB LIBSOURCES src/calc.cpp )
file(GLOB TESTS tests/calc_tests.cpp)

add_executable(CalcTests "test_runner.cpp" ${TESTS} ${LIBSOURCES})
target_link_libraries(CalcTests Catch)
target_include_directories(CalcTests PRIVATE include)