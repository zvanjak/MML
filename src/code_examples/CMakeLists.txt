project(MML_CodeExamples VERSION 1.0)

# Source files for concise code examples used by the documentation.
set(CODE_EXAMPLES_SOURCES
    code_examples_main.cpp
    readme00_fundamental_theorems.cpp
    readme01_quick_start.cpp
    readme02_vectors_matrices.cpp
    readme03_linear_systems.cpp
    readme04_functions_interpolation.cpp
    readme05_numerical_calculus.cpp
    readme06_field_operations.cpp
    readme07_ode_solvers.cpp
    readme09_polynomials.cpp
    readme10_root_finding.cpp
    readme11_coord_transforms.cpp
    readme12_parametric_curves.cpp
    readme13_path_integrals.cpp
    readme14_function_analysis.cpp
    readme15_dynamical_systems.cpp
)

add_executable(MML_CodeExamples ${CODE_EXAMPLES_SOURCES})

target_include_directories(MML_CodeExamples PRIVATE ${CMAKE_SOURCE_DIR})
