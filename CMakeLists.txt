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
    ${CMAKE_CURRENT_SOURCE_DIR}/src/book
    ${CMAKE_CURRENT_SOURCE_DIR}/src/book/mpl
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
add_subdirectory(src/book)
add_subdirectory(src/visualization_examples)
add_subdirectory(src/tryouts/test_linux_viz)
add_subdirectory(tools/utils_app)
add_subdirectory(src/examples)

# Book chapters (standalone applications)
add_subdirectory(src/book/chapters)

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
        test_linux_viz_standalone
        test_backend_switching
        test_all_visualizers
)

# If individual example executables exist, include them in the aggregate target
foreach(example_target
    Example00_ShockWave
    Example01_ProjectileLaunch
    Example02_DoublePendulum
    Example03_CollisionSimulator2D
    Example04_NBodyGravity
    Example05_LorentzTransformations
)
    if(TARGET ${example_target})
        add_dependencies(examples ${example_target})
    endif()
endforeach()

