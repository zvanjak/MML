cmake_minimum_required(VERSION 3.20)

project(MinimalMathLibrary VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Enable multi-processor compilation for MSVC (parallel within single target)
if(MSVC)
    add_compile_options(/MP)
endif()

# Options
option(MML_BUILD_DEVELOPMENT_TARGETS "Build demos, examples, tools, book apps, and tests" ON)
set(_MML_DEVELOPMENT_TARGET_DEFAULT ${MML_BUILD_DEVELOPMENT_TARGETS})
option(BUILD_SINGLE_HEADER_TESTS "Build single header tests (requires regenerating MML.h)" OFF)
option(BUILD_INTEGRATION_TESTS "Build long-running integration tests (N-body, large systems)" OFF)
option(MML_INSTALL "Generate install target" ON)
option(MML_BUILD_TESTS "Build MML unit test targets" ${_MML_DEVELOPMENT_TARGET_DEFAULT})
option(MML_BUILD_ALL_PRECISION_TESTS "Build test targets for float, double, and long double" OFF)
option(MML_BUILD_EXAMPLES "Build MML demo, README, and example applications" ${_MML_DEVELOPMENT_TARGET_DEFAULT})
option(MML_BUILD_DOCS_DEMOS "Build documentation demo applications" ${_MML_DEVELOPMENT_TARGET_DEFAULT})
option(MML_BUILD_BOOK "Build book and MPL demo applications" OFF)
option(MML_BUILD_VISUALIZATION_EXAMPLES "Build visualization example applications" ${_MML_DEVELOPMENT_TARGET_DEFAULT})
option(MML_BUILD_TOOLS "Build utility applications" OFF)
option(BUILD_TESTING "Build tests and enable CTest integration" ${_MML_DEVELOPMENT_TARGET_DEFAULT})

# =============================================================================
# Real type selection: float (default), double, or long_double
# Usage: cmake -B build -DMML_REAL_TYPE=float
# =============================================================================
set(MML_REAL_TYPE "float" CACHE STRING "Floating-point type for Real (float, double, long_double)")
set_property(CACHE MML_REAL_TYPE PROPERTY STRINGS "float" "double" "long_double")

if(MML_REAL_TYPE STREQUAL "float")
    add_compile_definitions(MML_USE_FLOAT)
    message(STATUS "MML Real type: float")
elseif(MML_REAL_TYPE STREQUAL "long_double")
    add_compile_definitions(MML_USE_LONG_DOUBLE)
    message(STATUS "MML Real type: long double")
elseif(MML_REAL_TYPE STREQUAL "double")
    add_compile_definitions(MML_USE_DOUBLE)
    message(STATUS "MML Real type: double")
else()
    message(FATAL_ERROR "Invalid MML_REAL_TYPE '${MML_REAL_TYPE}'. Use: float, double, or long_double")
endif()

# =============================================================================
# MML - Header-only INTERFACE library for package manager support
# =============================================================================
add_library(minimalmathlib INTERFACE)
add_library(minimalmathlib::minimalmathlib ALIAS minimalmathlib)

target_include_directories(minimalmathlib INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/mml/single_header>
    $<INSTALL_INTERFACE:include>

)

target_compile_features(minimalmathlib INTERFACE cxx_std_20)

# Global include directories for internal builds
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/mml/single_header
)

set(_MML_TESTS_ENABLED OFF)
if(BUILD_TESTING AND MML_BUILD_TESTS)
    set(_MML_TESTS_ENABLED ON)
endif()

set(_MML_INTEGRATION_TESTS_ENABLED OFF)
if(BUILD_TESTING AND BUILD_INTEGRATION_TESTS)
    set(_MML_INTEGRATION_TESTS_ENABLED ON)
endif()

if(_MML_TESTS_ENABLED OR _MML_INTEGRATION_TESTS_ENABLED)
    # Enable testing and generate CTest dashboard configuration.
    include(CTest)

    # Fetch Catch2 v3 for testing
    include(FetchContent)
    FetchContent_Declare(
      Catch2
      GIT_REPOSITORY https://github.com/catchorg/Catch2.git
      GIT_TAG        v3.11.0
    )
    FetchContent_MakeAvailable(Catch2)
    list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
endif()

# Treat diagnostics from first-party targets as build failures. This is placed
# after fetched dependencies so their warning policies remain independent.
if(MSVC)
    add_compile_options(/WX)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-Werror)
endif()

# Add subdirectories for optional applications and tooling
if(MML_BUILD_EXAMPLES)
    add_subdirectory(src/code_examples)
    add_subdirectory(src/examples)
endif()

if(MML_BUILD_DOCS_DEMOS)
    add_subdirectory(src/docs_demos)
endif()

if(MML_BUILD_BOOK AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/src/book/CMakeLists.txt")
    add_subdirectory(src/book)
endif()

if(MML_BUILD_BOOK AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/src/book/mpl_demo_app/CMakeLists.txt")
    add_subdirectory(src/book/mpl_demo_app)
endif()

if(MML_BUILD_VISUALIZATION_EXAMPLES)
    add_subdirectory(src/visualization_examples)
endif()

if(MML_BUILD_TOOLS AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/tools/utils_app/CMakeLists.txt")
    add_subdirectory(tools/utils_app)
endif()

# Book chapters (standalone applications)
if(MML_BUILD_BOOK AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/src/book/chapters/CMakeLists.txt")
    add_subdirectory(src/book/chapters)
endif()

# Precision testing application
if(_MML_TESTS_ENABLED)
    add_subdirectory(src/testing_precision)
endif()

# Tests
if(_MML_TESTS_ENABLED)
    add_subdirectory(tests)
endif()

# Integration Tests (optional - long-running validation tests)
if(_MML_INTEGRATION_TESTS_ENABLED)
    add_subdirectory(integration_tests)
endif()

# Aggregate target to build all example-like executables at once
set(_MML_EXAMPLE_AGGREGATE_TARGETS)
foreach(example_aggregate_target
    MML_CodeExamples
    MML_DocsApp
    MPL_DemoApp
    MML_VisualizationApp
)
    if(TARGET ${example_aggregate_target})
        list(APPEND _MML_EXAMPLE_AGGREGATE_TARGETS ${example_aggregate_target})
    endif()
endforeach()

# If individual example executables exist, include them in the aggregate target
foreach(example_target
    Example00_ShockWave
    Example01_ProjectileLaunch
    Example02_DoublePendulum
    Example03_CollisionSimulator2D
    Example04_NBodyGravity
    Example06_LorentzTransform
    Example07_EMFields
    Example08_GRExamples
    Example09_TypedFormsGradient
    Example10_TypedFormsFlux
    Example11_VectorSpaces
)
    if(TARGET ${example_target})
        list(APPEND _MML_EXAMPLE_AGGREGATE_TARGETS ${example_target})
    endif()
endforeach()

if(_MML_EXAMPLE_AGGREGATE_TARGETS)
    add_custom_target(examples
        DEPENDS ${_MML_EXAMPLE_AGGREGATE_TARGETS}
    )
endif()

# =============================================================================
# Installation rules for package managers (vcpkg, Conan, system install)
# =============================================================================
if(MML_INSTALL)
    include(GNUInstallDirs)
    include(CMakePackageConfigHelpers)

    # Install the INTERFACE target
    install(TARGETS minimalmathlib
        EXPORT minimalmathlib-targets
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    )

    # Install header files
    install(DIRECTORY mml/
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/mml
        FILES_MATCHING PATTERN "*.h"
    )

    # Install single-header version
    install(FILES mml/single_header/MML.h
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    )

    # Install license
    install(FILES LICENSE.md
        DESTINATION ${CMAKE_INSTALL_DOCDIR}
    )

    # Export targets
    install(EXPORT minimalmathlib-targets
        FILE minimalmathlib-targets.cmake
        NAMESPACE minimalmathlib::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/minimalmathlib
    )

    # Generate and install package config
    configure_package_config_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/cmake/minimalmathlib-config.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/minimalmathlib-config.cmake
        INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/minimalmathlib
    )

    # Generate version file
    write_basic_package_version_file(
        ${CMAKE_CURRENT_BINARY_DIR}/minimalmathlib-config-version.cmake
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY SameMajorVersion
    )

    # Install config files
    install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/minimalmathlib-config.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/minimalmathlib-config-version.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/minimalmathlib
    )
endif()
