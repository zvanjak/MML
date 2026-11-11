project(MML_Tests VERSION 1.0)

# Catch2 v3 is fetched in root CMakeLists.txt via FetchContent

# Compile foundational headers in separate translation units, without Catch2 or
# an umbrella header supplying dependencies before the header under test.
set(MML_ISOLATED_HEADERS
    mml/MMLBase.h
    mml/MMLTypeDefs.h
    mml/base/PrintFormatting.h
    mml/base/ExactInteger.h
    mml/base/BaseUtils/AngleCoordUtils.h
    mml/core/Integration/IntegrationBase.h
    mml/core/Integration/Integration1D.h
    mml/core/Integration/GaussKronrod.h
    mml/core/Integration/GaussianQuadrature.h
    mml/core/Integration/Integration2D.h
    mml/core/Integration/Integration3D.h
    mml/core/Integration/Integration2DAdaptive.h
    mml/core/Integration/Integration3DAdaptive.h
    mml/core/Integration/IntegrationImproper.h
    mml/core/Integration/SurfaceIntegration.h
    mml/core/Integration/PathIntegration.h
    mml/core/Integration/IntegrationParallel.h
    mml/core/Integration/MonteCarloIntegration.h
    mml/core/Integration/MonteCarloParallel.h
    mml/algorithms/MatrixAlg/MatrixInversion.h
    mml/algorithms/MatrixSingularity.h
    mml/algorithms/Eigen/HessenbergReduction.h
    mml/algorithms/Eigen/EigenSolver.h
    mml/algorithms/Eigen/detail/SchurEigenvectors.h
    mml/core/DifferentialGeometry/Metric.h
    mml/algorithms/ImplicitNewtonMatrix.h
    mml/algorithms/ODESolvers/ODEStepperInfrastructure.h
    mml/MMLExceptions.h
    mml/algorithms/Fourier/Fourier.h
    mml/algorithms/Optimization/LP/LinearProgram.h
    mml/algorithms/Optimization/LP/SimplexTableau.h
    mml/algorithms/Optimization/LP/SimplexSolver.h
    mml/algorithms/Optimization/Constraints/BoundConstraints.h
    mml/algorithms/Optimization/Constraints/BoxConstrainedNelderMead.h
    mml/algorithms/Optimization/Constraints/AugmentedLagrangian.h
    mml/algorithms/Optimization/LinearProgramming.h
    mml/algorithms/Statistics.h
    mml/algorithms/Statistics/Histogram.h
    mml/algorithms/Statistics/Distributions.h
    mml/algorithms/Statistics/DiscreteDistributions.h
    mml/algorithms/CurveFitting.h
    mml/MMLVisualizators.h
    mml/base/Intervals.h
    mml/base/Vector/Vector.h
    mml/base/Vector/VectorN.h
    mml/base/Vector/VectorTypes2D.h
    mml/base/Vector/VectorTypes3D.h
    mml/base/Vector/VectorTypes4D.h
    mml/base/Matrix/Matrix.h
    mml/base/Matrix/MatrixTriDiag.h
    mml/base/Geometry/Geometry.h
    mml/base/Geometry/GeometryBase/GeometryPoints.h
    mml/base/Geometry/GeometryBase/GeometryArithmetic.h
    mml/base/Geometry/Geometry3D/Geometry3DPlane.h
    mml/base/Geometry/Geometry3D/Geometry3DSurfaces.h
    mml/core/Derivation/DerivationRealFunction.h
    mml/core/Derivation/DerivationComplex.h
    mml/core/ComplexAnalysis.h
    mml/core/SparseSolvers/IterativeSolverBase.h
    mml/core/SparseSolvers/ConjugateGradient.h
    mml/core/SparseSolvers/BiCGSTAB.h
    mml/core/SparseSolvers/GMRES.h
    mml/core/SparseSolvers/Preconditioners.h
    mml/core/SparseSolvers/IterativeSolvers.h
    mml/core/LinAlgEqSolvers/LinAlgDirect.h
    mml/core/LinAlgEqSolvers/LinAlgSVD.h
    mml/core/Derivation/DerivationParametricCurve.h
    mml/core/Derivation/DerivationScalarFunction.h
    mml/base/Geometry/GeometryBase/Geometry2DShapes.h
    mml/base/Geometry/GeometryBase/Geometry3DShapes.h
    mml/interfaces/IFunction.h
    mml/base/StandardFunctions.h
    mml/core/DifferentialGeometry/Chart.h
    mml/core/DifferentialGeometry/FieldOperations.h
    mml/core/DifferentialGeometry/FormIntegration.h
    mml/algorithms/Analyzers/FunctionsAnalyzer.h
    mml/algorithms/Analyzers/MatrixAnalyzer.h
    mml/tools/analyzers/FunctionAnalysisReports.h
    mml/tools/analyzers/MatrixAnalysisReports.h
    mml/tools/visualizers/VisualizerConfig.h
    mml/tools/visualizers/VisualizerFiles.h
    mml/tools/visualizers/VisualizerProcess.h
    mml/tools/visualizers/VisualizerFacade.h
    mml/tools/persistence/VectorIO.h
    mml/tools/Visualizer.h
    test_data/parametric_curves_defs.h
    test_data/parametric_surfaces_defs.h
    test_beds/parametric_curves_test_bed.h
    test_beds/parametric_surfaces_test_bed.h
    test_beds/test_case_metadata.h
    test_beds/linear_alg_eq_systems_test_bed.h
)
if(BUILD_SINGLE_HEADER_TESTS)
    list(APPEND MML_ISOLATED_HEADERS mml/single_header/MML.h)
endif()
set(MML_ISOLATED_SOURCES)
foreach(header IN LISTS MML_ISOLATED_HEADERS)
    string(MAKE_C_IDENTIFIER "${header}" source_name)
    set(source "${CMAKE_CURRENT_BINARY_DIR}/header_checks/${source_name}.cpp")
    file(GENERATE OUTPUT "${source}" CONTENT "#include <${header}>\n")
    list(APPEND MML_ISOLATED_SOURCES "${source}")
endforeach()
add_library(MML_HeaderChecks OBJECT ${MML_ISOLATED_SOURCES}
    algorithms/analyzers/analyzer_boundary_checks.cpp
    base/geometry/geometry_deprecation_checks.cpp
    base/vector/vector_boundary_checks.cpp
    base/exception_boundary_checks.cpp
    base/matrix_algorithm_boundary_checks.cpp
    tools/visualizer_boundary_checks.cpp)
target_sources(MML_HeaderChecks PRIVATE base/namespace_collision_checks.cpp base/scalar_configuration_checks.cpp)
target_sources(MML_HeaderChecks PRIVATE base/namespace_macro_collision_checks.cpp)
target_sources(MML_HeaderChecks PRIVATE base/reserved_member_checks.cpp)
target_sources(MML_HeaderChecks PRIVATE base/algebra/algebra_boundary_checks.cpp)
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set_source_files_properties(base/reserved_member_checks.cpp PROPERTIES
        COMPILE_OPTIONS "-Werror=reserved-identifier")
endif()
target_sources(MML_HeaderChecks PRIVATE core/showcase_boundary_checks.cpp)
target_link_libraries(MML_HeaderChecks PRIVATE minimalmathlib::minimalmathlib)
set_target_properties(MML_HeaderChecks PROPERTIES UNITY_BUILD OFF DISABLE_PRECOMPILE_HEADERS ON)
if(MSVC)
    set_source_files_properties(base/geometry/geometry_deprecation_checks.cpp PROPERTIES
        COMPILE_OPTIONS "/w14996;/we4996")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    set_source_files_properties(base/geometry/geometry_deprecation_checks.cpp PROPERTIES
        COMPILE_OPTIONS "-Werror=deprecated-declarations")
endif()

# =============================================================================
# MML Core Tests - Tests for the core MML library
# =============================================================================

###############################################################################
# BASE TESTS - Fundamental data structures and utilities
###############################################################################
file(GLOB BASE_TESTS 
    # Algebra foundations
    "base/algebra/field_matrix_algorithms_tests.cpp"
    "base/algebra/algebra_traits_tests.cpp"
    "base/algebra/finite_groups_tests.cpp"
    "base/algebra/group_action_tests.cpp"
    "base/algebra/lie_groups_tests.cpp"
    "base/algebra/lie_group_contract_tests.cpp"
    "base/algebra/modular_arithmetic_tests.cpp"
    "base/algebra/polynomial_finite_field_tests.cpp"
    "base/algebra/representation_tests.cpp"
    "base/algebra/rigid_motion_groups_tests.cpp"
    "base/algebra/permutation_tests.cpp"
    "base/base_math_formula_fixes_tests.cpp"
    "base/header_hygiene_tests.cpp"
    "base/bugfix_epic_z5pg_tests.cpp"
    "base/bugfix_epic_i5u3_tests.cpp"
    "base/baseutils_tests.cpp"
    "base/baseutils_contract_tests.cpp"
    "base/coordinate_scale_contract_tests.cpp"
    "base/quaternion_scale_contract_tests.cpp"
    "base/matrix_scale_contract_tests.cpp"
    "base/interpolated_functions/construction_contract_tests.cpp"
    "base/interpolated_functions/numerical_contract_tests.cpp"
    "base/complex_instantiation_tests.cpp"
    "base/degenerate_inputs_tests.cpp"
    "base/matrix_naming_tests.cpp"
    "base/vector_naming_tests.cpp"
    "base/result_policy_tests.cpp"
    "base/uncomputed_result_contract_tests.cpp"
    "base/tolerance_policy_tests.cpp"
    "base/tolerance_base_tests.cpp"
    "base/matrixops_tests.cpp"
    "base/chebyshev_approximation_tests.cpp"
    "base/dirac_delta_function_tests.cpp"
    "base/exception_tests.cpp"
    "base/reserved_member_behavior_tests.cpp"
    "base/function_tests.cpp"
    "base/intervals_tests.cpp"
    "base/interval_set_contract_tests.cpp"
    "base/combinatorics_tests.cpp"
    "base/exact_integer_regression_tests.cpp"
    "base/mml_concepts_tests.cpp"
    "base/number_theory_tests.cpp"
    "base/ode_system_tests.cpp"
    "base/trajectory_storage_tests.cpp"
    "base/storage_boundary_regressions_tests.cpp"
    "base/polynom_tests.cpp"
    "base/quaternion_tests.cpp"
    "base/random_tests.cpp"
    "base/rational_tests.cpp"
    "base/quasi_random_tests.cpp"
    "base/special_functions_tests.cpp"
    "base/special_function_tail_contract_tests.cpp"
    "base/standard_functions_tests.cpp"
    "base/integer_boundary_contract_tests.cpp"
    "base/complex_reciprocal_contract_tests.cpp"
    "base/standard_fallback_contract_tests.cpp"
    # Differential forms
    "base/differential_forms/differential_form_tests.cpp"
    "base/differential_forms/hodge_tests.cpp"
    "base/differential_forms/hodge_scaling_tests.cpp"
    # Geometry tests
    "base/geometry/bounding_volumes_tests.cpp"
    "base/geometry/geometry_tests.cpp"
    "base/geometry/geometry_2d_tests.cpp"
    "base/geometry/geometry_3d_tests.cpp"
    "base/geometry/shape_formula_mesh_regressions_tests.cpp"
    "base/geometry/plane_rectangle_scale_contract_tests.cpp"
    "base/geometry/plane_query_range_contract_tests.cpp"
    "base/geometry/line_segment_contract_tests.cpp"
    "base/geometry/triangle_contract_tests.cpp"
    "base/geometry/distance_contract_tests.cpp"
    "base/geometry/polygon_contract_tests.cpp"
    "base/geometry/boundary_body_contract_tests.cpp"
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
    "base/tensor/tensor_contract_regressions_tests.cpp"
    # Vector tests
    "base/vector/vector_tests.cpp"
    "base/vector/vector_resize_contract_tests.cpp"
    "base/vector/unsigned_vector_tests.cpp"
    "base/vector/vectorn_tests.cpp"
    "base/vector/vectorn_subtypes_tests.cpp"
    "base/vector/vector_types_tests.cpp"
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
    "core/coord_transf/coordinate_transform_contract_tests.cpp"
    "interfaces/icoord_transf_tests.cpp"
    "interfaces/iparametrized_tests.cpp"
    "interfaces/idiscrete_map_tests.cpp"
    "interfaces/interface_naming_tests.cpp"
    "interfaces/ijacobian_contract_tests.cpp"
    "interfaces/default_implementation_tests.cpp"
    "interfaces/interface_validation_contract_tests.cpp"
    "core/coord_transf/coord_transf_quaternion_tests.cpp"
    "core/coord_transf/coord_transf_tensor_tests.cpp"
    "core/coord_transf/coord_transf_vector_tests.cpp"
    # Curves and surfaces
    "core/curves_tests.cpp"
    "core/curve_surface_contract_tests.cpp"
    "core/surface_step_contract_tests.cpp"
    "core/surfaces_tests.cpp"
    "core/showcase_geometry_tests.cpp"
    "core/differential_geometry/atlas_tests.cpp"
    "core/differential_geometry/chart_tests.cpp"
    "core/differential_geometry/form_integration_tests.cpp"
    "core/differential_geometry/frames_tests.cpp"
    "core/differential_geometry/induced_metric_tests.cpp"
    # Derivation
    "core/derivation/ad_consistency_tests.cpp"
    "core/derivation/ad_jacobians_tests.cpp"
    "core/derivation/ad_adapters_tests.cpp"
    "core/derivation/derivation_naming_tests.cpp"
    "core/derivation/derivation_characterization_tests.cpp"
    "core/derivation/derivation_tests.cpp"
    "core/derivation/richardson_tensor_contract_tests.cpp"
    "core/derivation/higher_derivative_error_contract_tests.cpp"
    "core/derivation/second_derivative_step_contract_tests.cpp"
    "core/derivation/jacobian_contract_tests.cpp"
    "core/derivation/derivative_step_contract_tests.cpp"
    "core/derivation/ad_quotient_constant_contract_tests.cpp"
    "core/derivation/ad_domain_type_contract_tests.cpp"
    "core/derivation/ad_tape_lifetime_contract_tests.cpp"
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
    # Integration
    "core/integration/integration_tests.cpp"
    "core/integration/integration_contract_regressions_tests.cpp"
    "core/integration/geometric_integration_contract_tests.cpp"
    "core/integration/integration_2d_tests.cpp"
    "core/integration/integration_3d_tests.cpp"
    "core/integration/gauss_kronrod_tests.cpp"
    "core/integration/gauss_kronrod_reference_tests.cpp"
    "core/integration/adaptive_integration_nd_tests.cpp"
    "core/integration/monte_carlo_integration_tests.cpp"
    "core/integration/parallel_integration_tests.cpp"
    "core/integration/path_integration_tests.cpp"
    "core/integration/surface_integration_tests.cpp"
    # Linear algebra
    "core/lin_alg_solvers/iterative_linear_solvers_tests.cpp"
    "core/lin_alg_solvers/dense_iterative_contract_tests.cpp"
    "core/lin_alg_solvers/linear_alg_eq_solvers_tests.cpp"
    "core/lin_alg_solvers/test_bed_metadata_tests.cpp"
    "base/vector/vector_unsigned_size_tests.cpp"
    "core/lin_alg_solvers/solver_scale_regressions_tests.cpp"
     "core/sparse_solvers/sparse_solvers_tests.cpp"
     "core/sparse_solvers/krylov_solvers_tests.cpp"
     "core/sparse_solvers/solver_boundary_contract_tests.cpp"
     "core/result_honesty_tests.cpp"
    "core/tolerance_core_tests.cpp"
    "core/metric_tensor_tests.cpp"
    "core/metric_field_contract_tests.cpp"
    "core/singularity_policy_tests.cpp"
    "core/mathematical_identity_tests.cpp"
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
    "core/orthogonal_basis_contract_tests.cpp"
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
    "algorithms/least_squares_contract_tests.cpp"
    "algorithms/generalized_least_squares_tests.cpp"
    "algorithms/domain_exception_tests.cpp"
    "algorithms/parametric_curve_analyzer_tests.cpp"
    # Eigensolvers
    "algorithms/complex_eigenanalysis_tests.cpp"
    "algorithms/eigensolver_adversarial_tests.cpp"
    "algorithms/eigensolver_building_blocks_tests.cpp"
    "algorithms/eigensystem_solvers_tests.cpp"
    "algorithms/self_adjoint_config_contract_tests.cpp"
    # Field analysis
    "algorithms/analyzers/field_analyzers_tests.cpp"
    "algorithms/analyzers/field_analyzer_contract_tests.cpp"
    "algorithms/analyzers/function_analyzer_tests.cpp"
    "algorithms/analyzers/function_analyzer_macro_contract_tests.cpp"
    "algorithms/analyzers/matrix_analyzer_tests.cpp"
    "algorithms/naming_honesty_tests.cpp"
    "algorithms/tolerance_algorithms_tests.cpp"
    # Root finding
    "algorithms/root_finding/root_finding_tests.cpp"
    "algorithms/root_finding/advanced_root_finding_tests.cpp"
    "algorithms/root_finding/convergence_contract_tests.cpp"
    "algorithms/root_finding/scalar_sign_contract_tests.cpp"
    "algorithms/callback_accounting_contract_tests.cpp"
    "algorithms/root_finding/polynomial_roots_tests.cpp"
    "algorithms/root_finding/polynomial_scale_contract_tests.cpp"
    "algorithms/root_finding/nonlinear_system_solvers_tests.cpp"
    # ODE solvers
    "algorithms/ode_solvers/ode_system_solvers_tests.cpp"
    "algorithms/ode_solvers/adaptive_acceptance_contract_tests.cpp"
    "algorithms/ode_solvers/adaptive_output_contract_tests.cpp"
    "algorithms/ode_solvers/adaptive_work_contract_tests.cpp"
    "algorithms/ode_solvers/stepper_infrastructure_header_tests.cpp"
    "algorithms/ode_solvers/bulirsch_stoer_extrapolation_tests.cpp"
    "algorithms/ode_solvers/stiff_benchmark_tests.cpp"
    "algorithms/ode_solvers/implicit_retry_tests.cpp"
    "algorithms/ode_solvers/trajectory_resource_tests.cpp"
    "algorithms/ode_solvers/reverse_time_tests.cpp"
    "algorithms/ode_solvers/ode_step_calculators_tests.cpp"
    "algorithms/ode_solvers/ode_event_detection_tests.cpp"
    "algorithms/ode_solvers/bvp_shooting_method_tests.cpp"
    "algorithms/ode_solvers/boundary_newton_status_tests.cpp"
    # DAE solvers
    "algorithms/dae_solvers/dae_solver_tests.cpp"
    "algorithms/dae_solvers/adaptive_dae_contract_tests.cpp"
    "algorithms/dae_solvers/linearly_implicit_contract_tests.cpp"
    "algorithms/dae_solvers/bdf4_contract_tests.cpp"
    "algorithms/dae_solvers/rodas4_tests.cpp"
    "algorithms/dae_solvers/dae_event_detection_tests.cpp"
    # Matrix algorithms
    "algorithms/matrix_alg_tests.cpp"
    "algorithms/faddeev_leverrier_contract_tests.cpp"
    "algorithms/matrix_inversion_consolidation_tests.cpp"
    # Graph algorithms
    "algorithms/graph_tests.cpp"
    "algorithms/graph_contract_tests.cpp"
    # Computational geometry
    "algorithms/comp_geometry/circles_tests.cpp"
    "algorithms/comp_geometry/convex_hull_tests.cpp"
    "algorithms/comp_geometry/convex_hull_3d_tests.cpp"
    "algorithms/comp_geometry/intersections_tests.cpp"
    "algorithms/comp_geometry/kdtree_tests.cpp"
    "algorithms/comp_geometry/kdtree_contract_tests.cpp"
    "algorithms/comp_geometry/polygon_ops_tests.cpp"
    "algorithms/comp_geometry/triangulation_tests.cpp"
    "algorithms/comp_geometry/triangulation_contract_tests.cpp"
    "algorithms/comp_geometry/circle_voronoi_contract_tests.cpp"
    "algorithms/comp_geometry/robust_predicate_contract_tests.cpp"
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
    "algorithms/optimization/lp_duality_contract_tests.cpp"
    "algorithms/optimization/lp_model_contract_tests.cpp"
    "algorithms/optimization/lp_sensitivity_contract_tests.cpp"
    "algorithms/optimization/augmented_lagrangian_tests.cpp"
    "algorithms/optimization/alm_contract_tests.cpp"
    # Statistics tests
    "algorithms/statistics/statistics_tests.cpp"
    "algorithms/statistics/statistics_scaling_contract_tests.cpp"
    "algorithms/statistics/probability_tail_contract_tests.cpp"
    "algorithms/statistics/quantile_contract_tests.cpp"
    "algorithms/statistics/discrete_boundary_contract_tests.cpp"
    "algorithms/statistics/distribution_parameter_contract_tests.cpp"
    "algorithms/statistics/discrete_sampling_contract_tests.cpp"
    "algorithms/statistics/distributions_tests.cpp"
    "algorithms/statistics/distribution_sampling_tests.cpp"
    "algorithms/statistics/discrete_distributions_tests.cpp"
    "algorithms/statistics/histogram_tests.cpp"
    # Fourier tests
    "algorithms/fourier/fourier_tests.cpp"
    "algorithms/fourier/dct_reference_contract_tests.cpp"
    "algorithms/fourier/fourier_detailed_contract_tests.cpp"
    "algorithms/fourier/fourier_helper_contract_tests.cpp"
    "algorithms/fourier/spectrum_layering_tests.cpp"
)

###############################################################################
# SYSTEMS TESTS - Linear algebra facade and dynamical systems
###############################################################################
file(GLOB SYSTEMS_TESTS
    "systems/linear_system_tests.cpp"
    "systems/linear_dispatch_contract_tests.cpp"
    "systems/linear_rhs_contract_tests.cpp"
    "systems/linear_verification_regressions_tests.cpp"
    "systems/dynamical_system_tests.cpp"
    "systems/fixed_point_contract_tests.cpp"
    "systems/lyapunov_contract_tests.cpp"
    "systems/systems_evaluation_contract_tests.cpp"
    "systems/parameter_scan_contract_tests.cpp"
    "systems/bifurcation_clock_contract_tests.cpp"
    "systems/report_qualification_contract_tests.cpp"
    "systems/model_domain_contract_tests.cpp"
)

###############################################################################
# TOOLS TESTS - Utility tools
###############################################################################
file(GLOB TOOLS_TESTS
    "tools/console_printer_tests.cpp"
    "tools/column_format_contract_tests.cpp"
    "tools/table_export_contract_tests.cpp"
    "tools/vector_table_contract_tests.cpp"
    "tools/shared_print_formatting_tests.cpp"
    "tools/analyzer_reports_tests.cpp"
    "tools/dataloader_tests.cpp"
    "tools/dataloader_numeric_contract_tests.cpp"
    "tools/dataloader_csv_record_tests.cpp"
    "tools/dataloader_case_folding_contract_tests.cpp"
    "tools/serializer_tests.cpp"
    "tools/serializer_preflight_contract_tests.cpp"
    "tools/serializer_metadata_contract_tests.cpp"
    "tools/serializer_sampling_contract_tests.cpp"
    "tools/persistence_portability_tests.cpp"
    "tools/binary_schema_contract_tests.cpp"
    "tools/vector_persistence_tests.cpp"
    "tools/json_conversion_contract_tests.cpp"
    "tools/persistence_save_contract_tests.cpp"
    "tools/particle_history_contract_tests.cpp"
    "tools/collection_contract_tests.cpp"
    "tools/timer_tests.cpp"
    "tools/scoped_timer_contract_tests.cpp"
    "tools/threadpool_tests.cpp"
    "tools/visualizer_tests.cpp"
    "tools/visualizer_posix_contract_tests.cpp"
    "tools/visualization_example_results_tests.cpp"
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
add_executable(MML_SolverBoundaryConsumerChecks EXCLUDE_FROM_ALL core/sparse_solvers/solver_boundary_consumer_checks.cpp)
target_link_libraries(MML_SolverBoundaryConsumerChecks PRIVATE minimalmathlib::minimalmathlib)
set_target_properties(MML_SolverBoundaryConsumerChecks PROPERTIES UNITY_BUILD OFF DISABLE_PRECOMPILE_HEADERS ON)
add_test(NAME MML_SolverBoundaryConsumerChecks COMMAND MML_SolverBoundaryConsumerChecks)
add_dependencies(MML_Tests MML_HeaderChecks MML_SolverBoundaryConsumerChecks)

add_executable(MML_VisualizerProcessFixture tools/visualizer_process_fixture.cpp)
add_executable(MML_ScopedTimerProcessFixture tools/scoped_timer_process_fixture.cpp)
target_link_libraries(MML_ScopedTimerProcessFixture PRIVATE minimalmathlib::minimalmathlib)
function(configure_mml_timer_fixture target fixture)
    add_dependencies(${target} ${fixture})
    target_compile_definitions(${target} PRIVATE MML_SCOPED_TIMER_FIXTURE="$<TARGET_FILE:${fixture}>")
endfunction()
function(configure_mml_launcher_fixture target)
    add_dependencies(${target} MML_VisualizerProcessFixture)
    target_compile_definitions(${target} PRIVATE
        MML_VISUALIZER_PROCESS_FIXTURE="$<TARGET_FILE:MML_VisualizerProcessFixture>")
endfunction()
configure_mml_launcher_fixture(MML_Tests)
configure_mml_timer_fixture(MML_Tests MML_ScopedTimerProcessFixture)

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
add_executable(MML_StandardFallback_Tests base/standard_fallback_contract_tests.cpp)
target_link_libraries(MML_StandardFallback_Tests PRIVATE Catch2::Catch2WithMain)
target_compile_definitions(MML_StandardFallback_Tests PRIVATE MML_HAS_STD_SPECIAL_FUNCTIONS=0 MML_EXPECT_STANDARD_FALLBACK)
catch_discover_tests(MML_StandardFallback_Tests)

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
#   & .\build\tests\Release\MML_Tests_Float.exe
#   & .\build\tests\Release\MML_Tests_LongDbl.exe
#   & .\build\tests\Release\MML_Tests.exe
# =============================================================================
function(configure_mml_test_precision target selection scalar)
    target_compile_definitions(${target} PRIVATE ${selection})
    foreach(other IN ITEMS MML_USE_FLOAT MML_USE_DOUBLE MML_USE_LONG_DOUBLE)
        if(NOT other STREQUAL selection)
            if(MSVC)
                target_compile_options(${target} PRIVATE "/U${other}")
            else()
                target_compile_options(${target} PRIVATE "-U${other}")
            endif()
        endif()
    endforeach()
    set(check "${CMAKE_CURRENT_BINARY_DIR}/${target}_precision.cpp")
    file(GENERATE OUTPUT "${check}" CONTENT
        "#include <mml/MMLBase.h>\n#include <type_traits>\nstatic_assert(std::is_same_v<MML::Real, ${scalar}>);\n")
    target_sources(${target} PRIVATE "${check}")
endfunction()

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
    configure_mml_launcher_fixture(MML_Tests_Float)
    configure_mml_timer_fixture(MML_Tests_Float MML_ScopedTimerProcessFixture)
    if(UNIX AND NOT APPLE)
        target_link_libraries(MML_Tests_Float PRIVATE atomic)
    endif()
    target_include_directories(MML_Tests_Float PRIVATE
        ${CMAKE_SOURCE_DIR}/mml/single_header
        ${CMAKE_SOURCE_DIR}/test_data
        ${CMAKE_SOURCE_DIR}/test_beds
        ${CMAKE_SOURCE_DIR}/tests
    )
    configure_mml_test_precision(MML_Tests_Float MML_USE_FLOAT float)
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
    configure_mml_launcher_fixture(MML_Tests_LongDbl)
    configure_mml_timer_fixture(MML_Tests_LongDbl MML_ScopedTimerProcessFixture)
    if(UNIX AND NOT APPLE)
        target_link_libraries(MML_Tests_LongDbl PRIVATE atomic)
    endif()
    target_include_directories(MML_Tests_LongDbl PRIVATE
        ${CMAKE_SOURCE_DIR}/mml/single_header
        ${CMAKE_SOURCE_DIR}/test_data
        ${CMAKE_SOURCE_DIR}/test_beds
        ${CMAKE_SOURCE_DIR}/tests
    )
    configure_mml_test_precision(MML_Tests_LongDbl MML_USE_LONG_DOUBLE "long double")
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
    add_executable(MML_SingleHeader_ScopedTimerProcessFixture tools/scoped_timer_process_fixture.cpp)
    target_include_directories(MML_SingleHeader_ScopedTimerProcessFixture PRIVATE ${CMAKE_SOURCE_DIR}/mml/single_header)
    target_compile_definitions(MML_SingleHeader_ScopedTimerProcessFixture PRIVATE MML_SINGLE_HEADER_INCLUDE=<MML.h>)
    configure_mml_timer_fixture(MML_SingleHeader_Tests MML_SingleHeader_ScopedTimerProcessFixture)
    target_sources(MML_SingleHeader_Tests PRIVATE tools/scoped_timer_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE
        base/coordinate_scale_contract_tests.cpp
        base/quaternion_scale_contract_tests.cpp
        base/matrix_scale_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE core/complex_svd_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE algorithms/faddeev_leverrier_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE algorithms/analyzers/field_analyzer_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE algorithms/analyzers/function_analyzer_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE algorithms/analyzers/function_analyzer_macro_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE algorithms/comp_geometry/kdtree_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE algorithms/comp_geometry/triangulation_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE algorithms/comp_geometry/circle_voronoi_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE algorithms/comp_geometry/robust_predicate_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE algorithms/analyzers/field_analyzers_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE tools/visualizer_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE tools/visualizer_posix_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE algorithms/graph_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE tools/vector_table_contract_tests.cpp)
    configure_mml_launcher_fixture(MML_SingleHeader_Tests)
    target_sources(MML_SingleHeader_Tests PRIVATE single_header/test_showcase_fixtures.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE single_header/namespace_collision_single_header.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE single_header/algebra_boundary_checks.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE base/matrix_naming_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE core/tolerance_core_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE algorithms/statistics/discrete_distributions_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE algorithms/ode_solvers/implicit_retry_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE algorithms/ode_solvers/boundary_newton_status_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE single_header/trajectory_storage_single_header.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE single_header/first_five_regressions_single_header.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE single_header/integration_contract_single_header.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE core/integration/integration_contract_regressions_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE core/integration/monte_carlo_integration_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE core/integration/geometric_integration_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE core/differential_geometry/form_integration_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE base/tensor/tensor_contract_regressions_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE base/geometry/shape_formula_mesh_regressions_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE base/geometry/geometry_spherical_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE base/geometry/plane_rectangle_scale_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE base/geometry/plane_query_range_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE base/geometry/line_segment_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE base/geometry/triangle_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE base/geometry/distance_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE base/geometry/polygon_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE base/geometry/boundary_body_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE base/differential_forms/hodge_scaling_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE systems/linear_dispatch_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE systems/linear_rhs_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE systems/linear_verification_regressions_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE algorithms/ode_solvers/adaptive_acceptance_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE algorithms/ode_solvers/adaptive_output_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE algorithms/ode_solvers/adaptive_work_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE algorithms/ode_solvers/bulirsch_stoer_extrapolation_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE algorithms/complex_eigenanalysis_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE algorithms/eigensolver_building_blocks_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE algorithms/eigensolver_adversarial_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE algorithms/self_adjoint_config_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE algorithms/root_finding/convergence_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE algorithms/root_finding/scalar_sign_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE algorithms/callback_accounting_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE algorithms/optimization/optimization_bound_constraints_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE algorithms/root_finding/polynomial_scale_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE tools/dataloader_numeric_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE tools/serializer_preflight_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE tools/serializer_metadata_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE tools/serializer_sampling_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE tools/binary_schema_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE tools/dataloader_csv_record_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE tools/dataloader_case_folding_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE tools/column_format_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE tools/table_export_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE algorithms/dae_solvers/adaptive_dae_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE algorithms/dae_solvers/linearly_implicit_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE algorithms/dae_solvers/bdf4_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE core/curve_surface_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE core/surface_step_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE core/metric_field_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE core/orthogonal_basis_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE core/lin_alg_solvers/dense_iterative_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE core/lin_alg_solvers/solver_scale_regressions_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE core/lin_alg_solvers/test_bed_metadata_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE base/vector/vector_unsigned_size_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE base/vector/vector_resize_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE core/derivation/richardson_tensor_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE core/derivation/higher_derivative_error_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE core/derivation/second_derivative_step_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE core/derivation/jacobian_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE core/derivation/derivative_step_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE core/derivation/ad_quotient_constant_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE core/derivation/ad_domain_type_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE core/derivation/ad_tape_lifetime_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE core/derivation/ad_adapters_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE algorithms/statistics/statistics_scaling_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE algorithms/statistics/probability_tail_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE algorithms/statistics/quantile_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE algorithms/statistics/discrete_boundary_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE algorithms/statistics/distribution_parameter_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE algorithms/statistics/discrete_sampling_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE tools/json_conversion_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE tools/persistence_save_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE tools/vector_persistence_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE tools/particle_history_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE tools/collection_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE systems/fixed_point_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE systems/lyapunov_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE systems/systems_evaluation_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE systems/parameter_scan_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE systems/bifurcation_clock_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE systems/report_qualification_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE systems/model_domain_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE base/special_function_tail_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE base/integer_boundary_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE interfaces/default_implementation_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE interfaces/interface_validation_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE base/interval_set_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE base/complex_reciprocal_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE base/standard_fallback_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE base/baseutils_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE base/uncomputed_result_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE base/interpolated_functions/construction_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE base/interpolated_functions/numerical_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE base/algebra/lie_group_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE algorithms/fourier/dct_reference_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE algorithms/fourier/fourier_detailed_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE algorithms/fourier/fourier_helper_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE core/complex_analysis_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE core/sparse_solvers/solver_boundary_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE algorithms/optimization/lp_duality_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE algorithms/optimization/lp_model_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE algorithms/optimization/lp_sensitivity_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE algorithms/optimization/alm_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE algorithms/analyzers/matrix_analyzer_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE single_header/generalized_least_squares_single_header.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE algorithms/least_squares_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE core/coord_transf/coordinate_transform_contract_tests.cpp)
    target_sources(MML_SingleHeader_Tests PRIVATE single_header/parallel_integration_single_header.cpp)
    target_include_directories(MML_SingleHeader_Tests PRIVATE ${CMAKE_SOURCE_DIR})
    add_executable(MML_SingleHeader_StandardFallback_Tests base/standard_fallback_contract_tests.cpp)
    target_link_libraries(MML_SingleHeader_StandardFallback_Tests PRIVATE Catch2::Catch2WithMain)
    target_include_directories(MML_SingleHeader_StandardFallback_Tests PRIVATE ${CMAKE_SOURCE_DIR}/mml/single_header)
    target_compile_definitions(MML_SingleHeader_StandardFallback_Tests PRIVATE MML_SINGLE_HEADER_INCLUDE=<MML.h> MML_HAS_STD_SPECIAL_FUNCTIONS=0 MML_EXPECT_STANDARD_FALLBACK)
    catch_discover_tests(MML_SingleHeader_StandardFallback_Tests)
    add_mml_single_header_test(MML_SingleHeader_Float_Tests MML_float.h)
    target_sources(MML_SingleHeader_Float_Tests PRIVATE core/complex_svd_tests.cpp)
    target_sources(MML_SingleHeader_Float_Tests PRIVATE base/geometry/geometry_spherical_tests.cpp)
    target_sources(MML_SingleHeader_Float_Tests PRIVATE base/differential_forms/hodge_scaling_tests.cpp)
    target_sources(MML_SingleHeader_Float_Tests PRIVATE base/matrix_naming_tests.cpp)
    target_sources(MML_SingleHeader_Float_Tests PRIVATE core/tolerance_core_tests.cpp)
    target_sources(MML_SingleHeader_Float_Tests PRIVATE algorithms/statistics/discrete_distributions_tests.cpp)
    target_sources(MML_SingleHeader_Float_Tests PRIVATE algorithms/statistics/quantile_contract_tests.cpp)
    target_sources(MML_SingleHeader_Float_Tests PRIVATE algorithms/ode_solvers/implicit_retry_tests.cpp)
    target_sources(MML_SingleHeader_Float_Tests PRIVATE algorithms/dae_solvers/adaptive_dae_contract_tests.cpp)
    target_sources(MML_SingleHeader_Float_Tests PRIVATE algorithms/statistics/discrete_sampling_contract_tests.cpp)
    target_sources(MML_SingleHeader_Float_Tests PRIVATE core/integration/geometric_integration_contract_tests.cpp)
    target_sources(MML_SingleHeader_Float_Tests PRIVATE core/derivation/jacobian_contract_tests.cpp)
    target_sources(MML_SingleHeader_Float_Tests PRIVATE algorithms/root_finding/convergence_contract_tests.cpp)
    target_sources(MML_SingleHeader_Float_Tests PRIVATE systems/bifurcation_clock_contract_tests.cpp)
    target_sources(MML_SingleHeader_Float_Tests PRIVATE algorithms/ode_solvers/adaptive_acceptance_contract_tests.cpp)
    target_sources(MML_SingleHeader_Float_Tests PRIVATE algorithms/ode_solvers/adaptive_output_contract_tests.cpp)
    target_sources(MML_SingleHeader_Float_Tests PRIVATE algorithms/ode_solvers/adaptive_work_contract_tests.cpp)
    target_sources(MML_SingleHeader_Float_Tests PRIVATE systems/lyapunov_contract_tests.cpp)
    target_sources(MML_SingleHeader_Float_Tests PRIVATE core/integration/integration_contract_regressions_tests.cpp)
    target_sources(MML_SingleHeader_Float_Tests PRIVATE base/special_function_tail_contract_tests.cpp)
    target_sources(MML_SingleHeader_Float_Tests PRIVATE algorithms/statistics/probability_tail_contract_tests.cpp)
    target_sources(MML_SingleHeader_Float_Tests PRIVATE tools/vector_persistence_tests.cpp)
    target_sources(MML_SingleHeader_Float_Tests PRIVATE base/geometry/distance_contract_tests.cpp)
    target_sources(MML_SingleHeader_Float_Tests PRIVATE core/derivation/ad_adapters_tests.cpp)
    target_sources(MML_SingleHeader_Float_Tests PRIVATE base/geometry/shape_formula_mesh_regressions_tests.cpp)
    target_sources(MML_SingleHeader_Float_Tests PRIVATE algorithms/ode_solvers/bulirsch_stoer_extrapolation_tests.cpp)
    target_sources(MML_SingleHeader_Float_Tests PRIVATE base/geometry/triangle_contract_tests.cpp)
    target_sources(MML_SingleHeader_Float_Tests PRIVATE systems/linear_verification_regressions_tests.cpp)
    target_compile_definitions(MML_SingleHeader_Float_Tests PRIVATE MML_SINGLE_HEADER_EXPECT_FLOAT)
    add_mml_single_header_test(MML_SingleHeader_LongDouble_Tests MML_long_double.h)
    target_sources(MML_SingleHeader_LongDouble_Tests PRIVATE core/complex_svd_tests.cpp)
    target_sources(MML_SingleHeader_LongDouble_Tests PRIVATE base/geometry/geometry_spherical_tests.cpp)
    target_sources(MML_SingleHeader_LongDouble_Tests PRIVATE base/differential_forms/hodge_scaling_tests.cpp)
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
