project(MML_DemoApp VERSION 1.0)

file(GLOB DEMO_APP_BASE_SOURCES
    "base/demo_baseutils.cpp"
    "base/demo_geometry.cpp"
    "base/demo_geometry_2d.cpp"
    "base/demo_geometry_3d.cpp"
    "base/demo_matrix.cpp"
    "base/demo_matrix_other.cpp"
    "base/demo_matrixnm.cpp"
    "base/demo_tensors.cpp"
    "base/demo_vector.cpp"
    "base/demo_vectorn.cpp"
    "base/demo_polynom.cpp"
    "base/demo_function.cpp"
    "base/demo_dirac_function.cpp"
)

file(GLOB DEMO_APP_CORE_SOURCES
    "core/demo_coord_transf.cpp"
    "core/demo_coord_system.cpp"
    "core/demo_curves_surfaces.cpp"
    "core/demo_derivation.cpp"
    "core/demo_integration.cpp"
    "core/demo_interpolated_function.cpp"
    "core/demo_linear_alg_eq_solvers.cpp"
    "core/demo_serialization.cpp"
    "core/demo_metric_tensor.cpp"
    "core/demo_field_operations.cpp"
    "core/demo_coord_system_old.cpp"
    "core/demo_coord_transf_vector.cpp"
    "core/demo_fields.cpp"
    "core/demo_intervals.cpp"
)

file(GLOB DEMO_APP_ALGORITHMS_SOURCES
    "algorithms/demo_diff_geometry.cpp"
    "algorithms/demo_function_analyzer.cpp"
    "algorithms/demo_ode_system_solvers.cpp"
    "algorithms/demo_path_integration.cpp"
    "algorithms/demo_root_finding.cpp"
    "algorithms/demo_statistics.cpp"
    "algorithms/demo_surface_integration.cpp"
)

file(GLOB TESTING_SPEED_SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/../testing_speed/*.cpp"
)

# Exclude main.cpp and test_framework_verify.cpp - they have their own main() functions
file(GLOB TESTING_PRECISION_SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/../testing_precision/test_precision_*.cpp"
)

# Note: Visualization examples are now a separate app (MML_VisualizationApp)
# See src/visualization_examples/CMakeLists.txt

add_executable(MML_DemoApp 
    demo_app_main.cpp
    ${DEMO_APP_BASE_SOURCES}
    ${DEMO_APP_CORE_SOURCES}
    ${DEMO_APP_ALGORITHMS_SOURCES}
    ${TESTING_SPEED_SOURCES}
    ${TESTING_PRECISION_SOURCES}
)

target_include_directories(MML_DemoApp PRIVATE ${CMAKE_SOURCE_DIR}/mml)
