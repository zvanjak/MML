cmake_minimum_required(VERSION 3.20)

project(MinimalMathLibrary VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Enable multi-processor compilation for MSVC (parallel within single target)
if(MSVC)
    add_compile_options(/MP)
endif()

# Global include directories for both MML and MPL libraries
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/mml
)

# Enable testing at the root level
enable_testing()

# Fetch Catch2 v3 for testing
include(FetchContent)
FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG        v3.11.0
)
FetchContent_MakeAvailable(Catch2)
list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)

# Add subdirectories for each application
add_subdirectory(src/demo_app)
add_subdirectory(src/docs_demos)
add_subdirectory(src/readme_examples)
add_subdirectory(src/visualization_examples)
add_subdirectory(src/examples)

# Precision testing application
add_subdirectory(src/testing_precision)

# Tests
add_subdirectory(tests)

# Aggregate target to build all example-like executables at once
add_custom_target(examples
    DEPENDS
        MML_ReadmeApp
        MML_DocsApp
        MML_DemoApp
        MML_VisualizationApp
)


