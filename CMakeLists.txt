cmake_minimum_required(VERSION 3.20)

project(MinimalMathLibrary VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Enable multi-processor compilation for MSVC (parallel within single target)
if(MSVC)
    add_compile_options(/MP)
endif()

# Options
option(BUILD_SINGLE_HEADER_TESTS "Build single header tests (requires regenerating MML.h)" OFF)
option(MML_INSTALL "Generate install target" ON)

# =============================================================================
# MML - Header-only INTERFACE library for package manager support
# =============================================================================
add_library(minimalmathlib INTERFACE)
add_library(minimalmathlib::minimalmathlib ALIAS minimalmathlib)

target_include_directories(minimalmathlib INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/mml>
    $<INSTALL_INTERFACE:include>
    $<INSTALL_INTERFACE:include/mml>
)

target_compile_features(minimalmathlib INTERFACE cxx_std_17)

# Global include directories for MML library (for internal builds)
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

# If individual example executables exist, include them in the aggregate target
foreach(example_target
    Example00_NBodyGravity
    Example01_ProjectileLaunch
    Example02_DoublePendulum
    Example03_F1GForces
    Example04_CollisionSim2D
    Example05_RigidBody
    Example06_LorentzTransform
)
    if(TARGET ${example_target})
        add_dependencies(examples ${example_target})
    endif()
endforeach()

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

