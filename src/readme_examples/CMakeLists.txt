project(MML_ReadmeApp VERSION 1.0)

# Explicitly list source files instead of GLOB to ensure CMake detects changes
set(README_APP_SOURCES
    readme_examples_main.cpp
    readme0_demonstration.cpp
    readme1_vectors.cpp
    readme2_matrices.cpp
    readme3_tensors.cpp
    readme4_geometry.cpp
    readme5_polynoms.cpp
    readme6_functions.cpp
    readme7_derivation.cpp
    readme8_integration.cpp
    readme9_vector_field_operations.cpp
    readme10_curves_surfaces.cpp
    readme11_coord_system_transf.cpp
    readme12_lin_alg_solvers.cpp
    readme13_ode_solvers.cpp
    readme14_function_analyzer.cpp
    readme15_visualizators.cpp
)

add_executable(MML_ReadmeApp ${README_APP_SOURCES})

target_include_directories(MML_ReadmeApp PRIVATE ${CMAKE_SOURCE_DIR}/mml)
