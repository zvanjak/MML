project(MML_Tests VERSION 1.0)

# Catch2 v3 is fetched in root CMakeLists.txt via FetchContent

# =============================================================================
# MML Core Tests - Tests for the core MML library
# =============================================================================

###############################################################################
# BASE TESTS - Fundamental data structures and utilities
###############################################################################
file(GLOB BASE_TESTS 
    # Algebra foundations
    "base/algebra/algebra_traits_tests.cpp"
    "base/algebra/finite_groups_tests.cpp"
    "base/algebra/group_action_tests.cpp"
    "base/algebra/lie_groups_tests.cpp"
    "base/algebra/modular_arithmetic_tests.cpp"
    "base/algebra/polynomial_finite_field_tests.cpp"
    "base/algebra/representation_tests.cpp"
    "base/algebra/rigid_motion_groups_tests.cpp"
    "base/algebra/permutation_tests.cpp"
    "base/base_math_formula_fixes_tests.cpp"
    "base/bugfix_epic_z5pg_tests.cpp"
    "base/bugfix_epic_i5u3_tests.cpp"
    "base/baseutils_tests.cpp"
    "base/matrixops_tests.cpp"
    "base/chebyshev_approximation_tests.cpp"
    "base/dirac_delta_function_tests.cpp"
    "base/exception_tests.cpp"
    "base/function_tests.cpp"
    "base/intervals_tests.cpp"
    "base/combinatorics_tests.cpp"
    "base/mml_concepts_tests.cpp"
    "base/number_theory_tests.cpp"
    "base/ode_system_tests.cpp"
    "base/polynom_tests.cpp"
    "base/quaternion_tests.cpp"
    "base/random_tests.cpp"
    "base/rational_tests.cpp"
    "base/quasi_random_tests.cpp"
    "base/special_functions_tests.cpp"
    "base/standard_functions_tests.cpp"
    # Differential forms
    "base/differential_forms/differential_form_tests.cpp"
    "base/differential_forms/hodge_tests.cpp"
    # Geometry tests
    "base/geometry/bounding_volumes_tests.cpp"
    "base/geometry/geometry_tests.cpp"
    "base/geometry/geometry_2d_tests.cpp"
    "base/geometry/geometry_3d_tests.cpp"
    "base/geometry/geometry_spherical_tests.cpp"
    "base/geometry/metric_musical_tests.cpp"
    "base/geometry/polygon_tests.cpp"
    # Geometry 3D bodies
    "base/geometry_3d_bodies/cube3d_tests.cpp"
    "base/geometry_3d_bodies/cubewithtriangles3d_tests.cpp"
    "base/geometry_3d_bodies/cylinder3d_tests.cpp"
    "base/geometry_3d_bodies/pyramid3d_tests.cpp"
    "base/geometry_3d_bodies/pyramidequilateral3d_tests.cpp"
    "base/geometry_3d_bodies/sphere3d_tests.cpp"
    "base/geometry_3d_bodies/torus3d_tests.cpp"
    # Interpolated functions
    "base/interpolated_functions/interpolated_functions_tests.cpp"
    # Matrix tests
    "base/matrix/matrix_tests.cpp"
    "base/matrix/matrix3d_tests.cpp"
    "base/matrix/matrix_sym_tests.cpp"
    "base/matrix/matrix_tridiag_tests.cpp"
    "base/matrix/matrix_banddiag_tests.cpp"
    "base/matrix/matrixNM_tests.cpp"
    "base/matrix/matrix_io_tests.cpp"
    # Sparse matrix tests (moved from MML Packages Private)
    "base/sparse_matrix/sparse_matrix_tests.cpp"
    # Tensor tests
    "base/tensor/rank1_tensor_tests.cpp"
    "base/tensor/tensors_tests.cpp"
    # Vector tests
    "base/vector/vector_tests.cpp"
    "base/vector/vectorn_tests.cpp"
    "base/vector/vectorn_subtypes_tests.cpp"
    "base/vector/vector_types_tests.cpp"
    "base/vector/vector_io_tests.cpp"
)

###############################################################################
# CORE TESTS - Core mathematical operations
###############################################################################
file(GLOB CORE_TESTS
    # Algebra operations
    "core/algebra/algebra_geometry_integration_tests.cpp"
    "core/algebra/algebra_law_checks_tests.cpp"
    "core/algebra/finite_group_algorithms_tests.cpp"
    "core/algebra/group_action_algorithms_tests.cpp"
    "core/algebra/modular_arithmetic_algorithms_tests.cpp"
    "core/algebra/polynomial_algorithms_tests.cpp"
    "core/algebra/representation_algorithms_tests.cpp"
    # Coordinate systems and transformations
    "core/coordinate_map_variance_tests.cpp"
    "core/coord_transf/coord_transf_tests.cpp"
    "core/coord_transf/coord_transf_quaternion_tests.cpp"
    "core/coord_transf/coord_transf_tensor_tests.cpp"
    "core/coord_transf/coord_transf_vector_tests.cpp"
    # Curves and surfaces
    "core/curves_tests.cpp"
    "core/surfaces_tests.cpp"
    "core/differential_geometry/atlas_tests.cpp"
    "core/differential_geometry/chart_tests.cpp"
    "core/differential_geometry/form_integration_tests.cpp"
    "core/differential_geometry/frames_tests.cpp"
    "core/differential_geometry/induced_metric_tests.cpp"
    # Derivation
    "core/derivation/ad_consistency_tests.cpp"
    "core/derivation/ad_jacobians_tests.cpp"
    "core/derivation/derivation_characterization_tests.cpp"
    "core/derivation/derivation_tests.cpp"
     "core/derivation/first_derivative_stencil_tests.cpp"
    "core/derivation/derivation_parametric_curve_tests.cpp"
    "core/derivation/derivation_parametric_surface_tests.cpp"
    "core/derivation/derivation_scalar_function_tests.cpp"
    "core/derivation/derivation_tensor_field_jacobians_tests.cpp"
    "core/derivation/derivation_vector_function_tests.cpp"
    "core/derivation/forward_ad_tests.cpp"
    "core/derivation/reverse_ad_tests.cpp"
    # Fields
    "core/fields/field_operations_tests.cpp"
    "core/fields/fields_tests.cpp"
    # Function helpers
    "core/function_helpers_tests.cpp"
    # Function spaces
    "core/function_spaces/function_spaces_tests.cpp"
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
    "core/lin_alg_solvers/iterative_linear_solvers_tests.cpp"
    "core/lin_alg_solvers/linear_alg_eq_solvers_tests.cpp"
     "core/sparse_solvers/sparse_solvers_tests.cpp"
    "core/metric_tensor_tests.cpp"
    "core/typed_field_operations_tests.cpp"
    "core/typed_metric_adapters_tests.cpp"
    "core/hermitian_cholesky_tests.cpp"
    "core/svd_tests.cpp"
    "core/complex_svd_tests.cpp"
    "core/vector_spaces/vector_spaces_tests.cpp"
    # Numeric validation
    "core/numeric_validation_tests.cpp"
    # Algorithm contract infrastructure
    "core/algorithm_types_tests.cpp"
    # Orthogonal basis
    "core/orthogonal_basis_tests.cpp"
    # Richardson extrapolation
    "core/richardson_extrapolation_tests.cpp"
    # Singularity handling
    "core/singularity_handling_tests.cpp"
    # Geometry/coordinate edge cases
    "core/coord_transf/geometry_coord_edge_case_tests.cpp"
    # Complex analysis
    "core/complex_analysis_tests.cpp"
)

###############################################################################
# ALGORITHMS TESTS - Numerical algorithms
###############################################################################
file(GLOB ALGORITHMS_TESTS
    # Curve fitting and analysis
    "algorithms/curve_fitting_tests.cpp"
    "algorithms/parametric_curve_analyzer_tests.cpp"
    # Eigensolvers
    "algorithms/complex_eigenanalysis_tests.cpp"
    "algorithms/eigensolver_building_blocks_tests.cpp"
    "algorithms/eigensystem_solvers_tests.cpp"
    # Field analysis
    "algorithms/analyzers/field_analyzers_tests.cpp"
    "algorithms/analyzers/function_analyzer_tests.cpp"
    "algorithms/analyzers/matrix_analyzer_tests.cpp"
    # Root finding
    "algorithms/root_finding/root_finding_tests.cpp"
    "algorithms/root_finding/advanced_root_finding_tests.cpp"
    "algorithms/root_finding/polynomial_roots_tests.cpp"
    "algorithms/root_finding/nonlinear_system_solvers_tests.cpp"
    # ODE solvers
    "algorithms/ode_solvers/ode_system_solvers_tests.cpp"
    "algorithms/ode_solvers/ode_step_calculators_tests.cpp"
    "algorithms/ode_solvers/ode_event_detection_tests.cpp"
    "algorithms/ode_solvers/bvp_shooting_method_tests.cpp"
    # DAE solvers
    "algorithms/dae_solvers/dae_solver_tests.cpp"
    "algorithms/dae_solvers/dae_event_detection_tests.cpp"
    # Matrix algorithms
    "algorithms/matrix_alg_tests.cpp"
    # Graph algorithms
    "algorithms/graph_tests.cpp"
    # Computational geometry
    "algorithms/comp_geometry/circles_tests.cpp"
    "algorithms/comp_geometry/convex_hull_tests.cpp"
    "algorithms/comp_geometry/convex_hull_3d_tests.cpp"
    "algorithms/comp_geometry/intersections_tests.cpp"
    "algorithms/comp_geometry/kdtree_tests.cpp"
    "algorithms/comp_geometry/polygon_ops_tests.cpp"
    "algorithms/comp_geometry/triangulation_tests.cpp"
    "algorithms/comp_geometry/voronoi_tests.cpp"
    # Optimization algorithms (1D minimization)
    "algorithms/optimization/optimization_tests.cpp"
    # Optimization algorithms (multidimensional)
    "algorithms/optimization/optimization_multidim_tests.cpp"
    # Optimization algorithms (bound constraints)
    "algorithms/optimization/optimization_bound_constraints_tests.cpp"
    "algorithms/optimization/optimization_box_constrained_nelder_mead_tests.cpp"
    "algorithms/optimization/optimization_box_constrained_powell_tests.cpp"
    "algorithms/optimization/optimization_projected_gradient_tests.cpp"
    # Optimization algorithms (global stochastic)
    # Optimization algorithms (linear programming)
    "algorithms/optimization/optimization_linear_programming_tests.cpp"
    "algorithms/optimization/augmented_lagrangian_tests.cpp"
    # Statistics tests
    "algorithms/statistics/statistics_tests.cpp"
    "algorithms/statistics/distributions_tests.cpp"
    "algorithms/statistics/distribution_sampling_tests.cpp"
    "algorithms/statistics/discrete_distributions_tests.cpp"
    "algorithms/statistics/histogram_tests.cpp"
    # Fourier tests
    "algorithms/fourier/fourier_tests.cpp"
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

# Link libatomic on Linux for 128-bit atomics (needed by long double builds)
if(UNIX AND NOT APPLE)
    target_link_libraries(MML_Tests PRIVATE atomic)
endif()

target_include_directories(MML_Tests PRIVATE 
    ${CMAKE_SOURCE_DIR}/mml/single_header
    ${CMAKE_SOURCE_DIR}/test_data
    ${CMAKE_SOURCE_DIR}/test_beds
    ${CMAKE_SOURCE_DIR}/tests
)

# Add /FS flag for MSVC to allow parallel compilation with shared PDB
if(MSVC)
    target_compile_options(MML_Tests PRIVATE /FS)
endif()

# Use catch_discover_tests for better CTest integration
# Label allows filtering: ctest -L core -j8
include(Catch)
catch_discover_tests(MML_Tests 
    PROPERTIES LABELS "core"
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

# =============================================================================
# Multi-Precision Test Targets (float / long double)
# =============================================================================
# Build separate test executables for each Real type, each with its own
# compile definition overriding the configured default.
# Enable with: cmake -B build -DMML_BUILD_ALL_PRECISION_TESTS=ON
#
# Usage:
#   ctest -L float          # Run only float tests
#   ctest -L long_double    # Run only long double tests
#   ctest -L core           # Run only tests for the configured default
# =============================================================================
if(MML_BUILD_ALL_PRECISION_TESTS)
    # --- Float tests ---
    add_executable(MML_Tests_Float
        ${BASE_TESTS}
        ${CORE_TESTS}
        ${ALGORITHMS_TESTS}
        ${SYSTEMS_TESTS}
        ${TOOLS_TESTS}
        ${INFRASTRUCTURE_TESTS}
    )
    target_link_libraries(MML_Tests_Float PRIVATE Catch2::Catch2WithMain)
    if(UNIX AND NOT APPLE)
        target_link_libraries(MML_Tests_Float PRIVATE atomic)
    endif()
    target_include_directories(MML_Tests_Float PRIVATE
        ${CMAKE_SOURCE_DIR}/mml/single_header
        ${CMAKE_SOURCE_DIR}/test_data
        ${CMAKE_SOURCE_DIR}/test_beds
        ${CMAKE_SOURCE_DIR}/tests
    )
    target_compile_definitions(MML_Tests_Float PRIVATE MML_USE_FLOAT)
    if(MSVC)
        target_compile_options(MML_Tests_Float PRIVATE /FS)
    endif()
    catch_discover_tests(MML_Tests_Float
        PROPERTIES LABELS "float"
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )

    # --- Long Double tests ---
    add_executable(MML_Tests_LongDbl
        ${BASE_TESTS}
        ${CORE_TESTS}
        ${ALGORITHMS_TESTS}
        ${SYSTEMS_TESTS}
        ${TOOLS_TESTS}
        ${INFRASTRUCTURE_TESTS}
    )
    target_link_libraries(MML_Tests_LongDbl PRIVATE Catch2::Catch2WithMain)
    if(UNIX AND NOT APPLE)
        target_link_libraries(MML_Tests_LongDbl PRIVATE atomic)
    endif()
    target_include_directories(MML_Tests_LongDbl PRIVATE
        ${CMAKE_SOURCE_DIR}/mml/single_header
        ${CMAKE_SOURCE_DIR}/test_data
        ${CMAKE_SOURCE_DIR}/test_beds
        ${CMAKE_SOURCE_DIR}/tests
    )
    target_compile_definitions(MML_Tests_LongDbl PRIVATE MML_USE_LONG_DOUBLE)
    if(MSVC)
        target_compile_options(MML_Tests_LongDbl PRIVATE /FS)
    endif()
    catch_discover_tests(MML_Tests_LongDbl
        PROPERTIES LABELS "long_double"
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )

    message(STATUS "Multi-precision tests enabled: MML_Tests_Float, MML_Tests_LongDbl")
endif()

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
    function(add_mml_single_header_test target_name header_name)
        add_executable(${target_name} single_header/test_single_header.cpp)
        target_link_libraries(${target_name} PRIVATE Catch2::Catch2WithMain)
        target_include_directories(${target_name} PRIVATE ${CMAKE_SOURCE_DIR}/mml/single_header)
        target_compile_definitions(${target_name} PRIVATE MML_SINGLE_HEADER_INCLUDE=<${header_name}>)
        catch_discover_tests(${target_name})
    endfunction()

    add_mml_single_header_test(MML_SingleHeader_Tests MML.h)
    add_mml_single_header_test(MML_SingleHeader_Float_Tests MML_float.h)
    target_compile_definitions(MML_SingleHeader_Float_Tests PRIVATE MML_SINGLE_HEADER_EXPECT_FLOAT)
    add_mml_single_header_test(MML_SingleHeader_LongDouble_Tests MML_long_double.h)
    target_compile_definitions(MML_SingleHeader_LongDouble_Tests PRIVATE MML_SINGLE_HEADER_EXPECT_LONG_DOUBLE)
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
target_include_directories(GenerateSerializerRefFiles PRIVATE ${CMAKE_SOURCE_DIR}/mml/single_header)

