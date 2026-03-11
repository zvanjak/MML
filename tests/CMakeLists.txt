project(MML_Tests VERSION 1.0)

# Catch2 v3 is fetched in root CMakeLists.txt via FetchContent

# =============================================================================
# MML Core Tests
# =============================================================================

###############################################################################
# BASE TESTS - Fundamental data structures and utilities
###############################################################################
file(GLOB BASE_TESTS 
    "base/base_math_formula_fixes_tests.cpp"
    "base/baseutils_tests.cpp"
    "base/chebyshev_approximation_tests.cpp"
    "base/dirac_delta_function_tests.cpp"
    "base/function_tests.cpp"
    # Geometry tests
    "base/geometry/bounding_volumes_tests.cpp"
    "base/geometry/geometry_tests.cpp"
    "base/geometry/geometry_2d_tests.cpp"
    "base/geometry/geometry_3d_tests.cpp"
    "base/geometry/geometry_spherical_tests.cpp"
    "base/geometry/polygon_tests.cpp"
    # Geometry 3D bodies
    "base/geometry_3d_bodies/cube3d_tests.cpp"
    "base/geometry_3d_bodies/cubewithtriangles3d_tests.cpp"
    "base/geometry_3d_bodies/cylinder3d_tests.cpp"
    "base/geometry_3d_bodies/pyramid3d_tests.cpp"
    "base/geometry_3d_bodies/pyramidequilateral3d_tests.cpp"
    "base/geometry_3d_bodies/sphere3d_tests.cpp"
    "base/geometry_3d_bodies/torus3d_tests.cpp"
    "base/interpolated_functions_tests.cpp"
    "base/intervals_tests.cpp"
    # Matrix tests
    "base/matrix/matrix_tests.cpp"
    "base/matrix/matrix3d_tests.cpp"
    "base/matrix/matrix_sym_tests.cpp"
    "base/matrix/matrix_tridiag_tests.cpp"
    "base/matrix/matrix_banddiag_tests.cpp"
    "base/matrix/matrixNM_tests.cpp"
    "base/matrix_io_tests.cpp"
    "base/ode_system_tests.cpp"
    "base/polynom_tests.cpp"
    "base/quaternion_tests.cpp"
    "base/random_tests.cpp"
    "base/standard_functions_tests.cpp"
    "base/tensors_tests.cpp"
    # Vector tests
    "base/vector/vector_tests.cpp"
    "base/vector/vectorn_tests.cpp"
    "base/vector/vectorn_subtypes_tests.cpp"
    "base/vector/vector_types_tests.cpp"
    "base/vector_io_tests.cpp"
)

###############################################################################
# CORE TESTS - Core mathematical operations
###############################################################################
file(GLOB CORE_TESTS
    # Coordinate systems and transformations
    "core/coord_system_tests.cpp"
    "core/coord_transf/coord_transf_tests.cpp"
    "core/coord_transf/coord_transf_lorentz_tests.cpp"
    "core/coord_transf/coord_transf_quaternion_tests.cpp"
    "core/coord_transf/coord_transf_tensor_tests.cpp"
    "core/coord_transf/coord_transf_vector_tests.cpp"
    # Curves and surfaces
    "core/curves_tests.cpp"
    "core/surfaces_tests.cpp"
    # Derivation
    "core/derivation/derivation_tests.cpp"
    "core/derivation/derivation_parametric_curve_tests.cpp"
    "core/derivation/derivation_parametric_surface_tests.cpp"
    "core/derivation/derivation_scalar_function_tests.cpp"
    "core/derivation/derivation_tensor_field_jacobians_tests.cpp"
    "core/derivation/derivation_vector_function_tests.cpp"
    # Field operations
    "core/field_operations_tests.cpp"
    # Fields
    "core/fields_tests.cpp"
    # Function helpers
    "core/function_helpers_tests.cpp"
    # Integration
    "core/integration/integration_tests.cpp"
    "core/integration/integration_2d_tests.cpp"
    "core/integration/integration_3d_tests.cpp"
    "core/integration/gauss_kronrod_tests.cpp"
    "core/integration/adaptive_integration_nd_tests.cpp"
    "core/integration/monte_carlo_integration_tests.cpp"
    "core/integration/path_integration_tests.cpp"
    "core/integration/surface_integration_tests.cpp"
    # Linear algebra
    "core/iterative_linear_solvers_tests.cpp"
    "core/linear_alg_eq_solvers_tests.cpp"
    "core/matrixutils_tests.cpp"
    "core/metric_tensor_tests.cpp"
    "core/svd_tests.cpp"
    # Numeric validation
    "core/numeric_validation_tests.cpp"
    # Orthogonal basis
    "core/orthogonal_basis_tests.cpp"
    # Richardson extrapolation
    "core/richardson_extrapolation_tests.cpp"
    # Singularity handling
    "core/singularity_handling_tests.cpp"
    # Geometry/coordinate edge cases
    "core/geometry_coord_edge_case_tests.cpp"
)

###############################################################################
# ALGORITHMS TESTS - Numerical algorithms
###############################################################################
file(GLOB ALGORITHMS_TESTS
    # Curve fitting and analysis
    "algorithms/curve_fitting_tests.cpp"
    "algorithms/parametric_curve_analyzer_tests.cpp"
    # Eigensolvers
    "algorithms/eigensolver_building_blocks_tests.cpp"
    "algorithms/eigensystem_solvers_tests.cpp"
    # Field analysis
    "algorithms/field_analyzers_tests.cpp"
    "algorithms/function_analyzer_tests.cpp"
    # Root finding
    "algorithms/root_finding_tests.cpp"
    "algorithms/advanced_root_finding_tests.cpp"
    "algorithms/polynomial_roots_tests.cpp"
    # ODE solvers
    "algorithms/ode_solvers/ode_system_solvers_tests.cpp"
    "algorithms/ode_solvers/ode_step_calculators_tests.cpp"
    "algorithms/ode_solvers/ode_event_detection_tests.cpp"
    "algorithms/bvp_shooting_method_tests.cpp"
    # DAE solvers
    "algorithms/dae_solver_tests.cpp"
    # Matrix algorithms
    "algorithms/matrix_alg_tests.cpp"
    # Graph algorithms
    "algorithms/graph_tests.cpp"
    # Computational geometry (split into subdirectory)
    "algorithms/comp_geometry/circles_tests.cpp"
    "algorithms/comp_geometry/convex_hull_tests.cpp"
    "algorithms/comp_geometry/convex_hull_3d_tests.cpp"
    "algorithms/comp_geometry/intersections_tests.cpp"
    "algorithms/comp_geometry/kdtree_tests.cpp"
    "algorithms/comp_geometry/polygon_ops_tests.cpp"
    "algorithms/comp_geometry/triangulation_tests.cpp"
    "algorithms/comp_geometry/voronoi_tests.cpp"
    # Genetic algorithms (core algorithm, also used by optimization package)
    "algorithms/genetic_algorithm_tests.cpp"
    # Optimization algorithms (1D minimization)
    "algorithms/optimization_tests.cpp"
    # Statistics tests
    "algorithms/statistics_tests.cpp"
)

###############################################################################
# SYSTEMS TESTS - Linear algebra facade and dynamical systems
###############################################################################
file(GLOB SYSTEMS_TESTS
    "systems/linear_system_tests.cpp"
    "systems/dynamical_system_tests.cpp"
)

###############################################################################
# TOOLS TESTS - Utility tools
###############################################################################
file(GLOB TOOLS_TESTS
    "tools/console_printer_tests.cpp"
    "tools/dataloader_tests.cpp"
    "tools/serializer_tests.cpp"
    "tools/timer_tests.cpp"
    "tools/threadpool_tests.cpp"
    "tools/visualizer_tests.cpp"
)

###############################################################################
# INFRASTRUCTURE TESTS - Test framework utilities
###############################################################################
set(INFRASTRUCTURE_TESTS
    "test_precision_helpers_tests.cpp"
    "test_matchers_tests.cpp"
)

###############################################################################
# BUILD MML_Tests - Core library tests only
###############################################################################
add_executable(MML_Tests 
    ${BASE_TESTS}
    ${CORE_TESTS}
    ${ALGORITHMS_TESTS}
    ${SYSTEMS_TESTS}
    ${TOOLS_TESTS}
    ${INFRASTRUCTURE_TESTS}
)

target_link_libraries(MML_Tests PRIVATE 
    Catch2::Catch2WithMain
)

target_include_directories(MML_Tests PRIVATE 
    ${CMAKE_SOURCE_DIR}/mml
    ${CMAKE_SOURCE_DIR}/test_data
    ${CMAKE_SOURCE_DIR}/test_beds
    ${CMAKE_SOURCE_DIR}/tests
)

# Add /FS flag for MSVC to allow parallel compilation with shared PDB
if(MSVC)
    target_compile_options(MML_Tests PRIVATE /FS)
endif()

# Precompiled header - dramatically reduces compile time by pre-parsing
# catch2 headers once instead of 100+ times
target_precompile_headers(MML_Tests PRIVATE
    <catch2/catch_all.hpp>
)

# Register test executable with CTest directly (simple, no discovery overhead)
add_test(NAME MML_Tests 
    COMMAND MML_Tests --reporter compact
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
set_tests_properties(MML_Tests PROPERTIES LABELS "core" TIMEOUT 600)

# ============================================================================
# Single Header Test - DISABLED pending generator fixes
# ============================================================================
# Issues found in MML.h generator (MinimalMathLibrary-m8rl):
# Single header tests - disabled by default because they require regenerating MML.h
# Enable with: cmake -B build -DBUILD_SINGLE_HEADER_TESTS=ON
#
# Known issues when single header generation is broken:
# 1. File ordering: MMLPrecision.h must come before MMLBase.h
# 2. BOM characters in some source files not stripped
# 3. Missing <filesystem> header
# 4. Class redefinitions (ScalarFieldAnalyzer, VectorFieldAnalyzer)
# 5. Forward reference issues (DCT, EigenSolver, etc.)
if(BUILD_SINGLE_HEADER_TESTS)
    add_executable(MML_SingleHeader_Tests single_header/test_single_header.cpp)
    target_link_libraries(MML_SingleHeader_Tests PRIVATE Catch2::Catch2WithMain)
    target_include_directories(MML_SingleHeader_Tests PRIVATE ${CMAKE_SOURCE_DIR}/mml)
    catch_discover_tests(MML_SingleHeader_Tests)
endif()

# ============================================================================
# Serializer Reference File Generator - Run once to generate golden files
# ============================================================================
# This utility generates reference files for serializer golden file tests.
# Run it manually when serialization format changes intentionally.
# Output: tests/tools/serialization_test_data/*.txt
add_executable(GenerateSerializerRefFiles 
    tools/generate_serializer_reference_files.cpp
)
target_include_directories(GenerateSerializerRefFiles PRIVATE ${CMAKE_SOURCE_DIR}/mml)

