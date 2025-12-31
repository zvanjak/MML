project(MML_Tests VERSION 1.0)

# Catch2 v3 is fetched in root CMakeLists.txt via FetchContent

# Collect all test files
# NOTE: Most tests commented out during Real type flexibility implementation
# Uncomment progressively as precision-aware matchers are implemented
file(GLOB BASE_TESTS 
    "base/baseutils_tests.cpp"          # ACTIVE - migrated to precision-aware
    "base/bounding_volumes_tests.cpp"    # ACTIVE - migrated to precision-aware
    "base/chebyshev_approximation_tests.cpp"  # ACTIVE - migrated to precision-aware
    "base/dirac_delta_function_tests.cpp"  # ACTIVE - migrated to precision-aware
    "base/geometry_tests.cpp"           # ACTIVE - migrated to precision-aware
    "base/geometry_2d_tests.cpp"        # ACTIVE - migrated to precision-aware
    "base/geometry_3d_tests.cpp"        # ACTIVE - migrated to precision-aware
    "base/geometry_spherical_tests.cpp"
    "base/geometry_3d_bodies/cube3d_tests.cpp"
    "base/geometry_3d_bodies/cubewithtriangles3d_tests.cpp"
    "base/geometry_3d_bodies/cylinder3d_tests.cpp"
    "base/geometry_3d_bodies/pyramid3d_tests.cpp"
    "base/geometry_3d_bodies/pyramidequilateral3d_tests.cpp"
    "base/geometry_3d_bodies/sphere3d_tests.cpp"
    "base/geometry_3d_bodies/torus3d_tests.cpp"
    "base/intervals_tests.cpp"          # ACTIVE - migrated to precision-aware
    "base/matrix_tests.cpp"             # ACTIVE - migrated to precision-aware
    "base/matrix_sym_tests.cpp"         # ACTIVE - migrated to precision-aware
    "base/matrix_tridiag_tests.cpp"     # ACTIVE - migrated to precision-aware
    "base/matrix_banddiag_tests.cpp"    # ACTIVE - migrated to precision-aware
    "base/matrixNM_tests.cpp"           # ACTIVE - migrated to precision-aware
    "base/ode_system_tests.cpp"         # ACTIVE - migrated to precision-aware
    "base/polygon_tests.cpp"            # ACTIVE - migrated to precision-aware
    "base/tensors_tests.cpp"            # ACTIVE - migrated to precision-aware
    "base/vector_tests.cpp"             # ACTIVE - basic vectors
    "base/vectorN_tests.cpp"            # ACTIVE - migrated to precision-aware
    "base/vectorN_subtypes_tests.cpp"   # ACTIVE - migrated to precision-aware
    "base/polynom_tests.cpp"            # ACTIVE - migrated to precision-aware
)

file(GLOB CORE_TESTS
    "core/coord_transf_tests.cpp"       # ACTIVE - migrated to precision-aware
    "core/coord_transf_quaternion_tests.cpp"  # ACTIVE - migrated to precision-aware
    "core/coord_transf_lorentz_tests.cpp"  # ACTIVE - migrated to precision-aware
    "core/coord_transf_vector_tests.cpp" # ACTIVE - migrated to precision-aware
    "core/curves_tests.cpp"             # ACTIVE - migrated to precision-aware
    "core/surfaces_tests.cpp"           # ACTIVE - migrated to precision-aware
    "core/derivation_tests.cpp"         # ACTIVE - migrated to precision-aware
    "core/derivation_parametric_curve_tests.cpp"  # ACTIVE - migrated to precision-aware
    "core/derivation_parametric_surface_tests.cpp"  # ACTIVE - migrated to precision-aware
    "core/derivation_scalar_function_tests.cpp"  # ACTIVE - migrated to precision-aware
    "core/derivation_tensor_field_jacobians_tests.cpp" # ACTIVE - migrated to precision-aware
    "core/derivation_vector_function_tests.cpp" # ACTIVE - migrated to precision-aware
    "core/field_operations_tests.cpp"       # ACTIVE - migrated to precision-aware
    "core/function_tests.cpp"           # ACTIVE - migrated to precision-aware
    "core/integration_tests.cpp"        # ACTIVE - migrated to precision-aware
    "core/interpolated_functions_tests.cpp" # ACTIVE - migrated to precision-aware
    "core/iterative_linear_solvers_tests.cpp" # ACTIVE - migrated to precision-aware
    "core/linear_alg_eq_solvers_tests.cpp"  # ACTIVE - migrated to precision-aware
    "core/matrixutils_tests.cpp"            # ACTIVE - migrated to precision-aware
    "core/metric_tensor_tests.cpp"          # ACTIVE - migrated to precision-aware
    "core/svd_tests.cpp"                    # ACTIVE - migrated to precision-aware
)

file(GLOB ALGORITHMS_TESTS
    "algorithms/curve_fitting_tests.cpp"       # ACTIVE - migrated to precision-aware
    "algorithms/eigensystem_solvers_tests.cpp" # ACTIVE - migrated to precision-aware
    "algorithms/field_analyzers_tests.cpp"      # ACTIVE - migrated to precision-aware
    "algorithms/function_analyzer_tests.cpp"    # ACTIVE - migrated to precision-aware
    "algorithms/integration_2d_tests.cpp"       # ACTIVE - 2D double integration tests
    "algorithms/matrix_alg_tests.cpp"           # ACTIVE - migrated to precision-aware
    "algorithms/monte_carlo_integration_tests.cpp" # ACTIVE - migrated to precision-aware
    "algorithms/ode_system_solvers_tests.cpp"   # ACTIVE - migrated to precision-aware
    "algorithms/optimization_tests.cpp"         # ACTIVE - migrated to precision-aware
    "algorithms/optimization_multidim_tests.cpp"  # ACTIVE - migrated to precision-aware
    "algorithms/optimization_heuristic_tests.cpp" # ACTIVE - migrated to precision-aware
    "algorithms/parametric_curve_analyzer_tests.cpp" # ACTIVE - migrated to precision-aware
    "algorithms/path_integration_tests.cpp"     # ACTIVE - migrated to precision-aware
    "algorithms/polynomial_roots_tests.cpp"     # ACTIVE - migrated to precision-aware
    "algorithms/quaternion_tests.cpp"           # ACTIVE - migrated to precision-aware
    "algorithms/root_finding_tests.cpp"         # ACTIVE - migrated to precision-aware
    "algorithms/advanced_root_finding_tests.cpp" # ACTIVE - migrated to precision-aware
    "algorithms/statistics_tests.cpp"            # ACTIVE - migrated to precision-aware
    "algorithms/time_series_tests.cpp"          # ACTIVE - time series analysis
    "algorithms/histogram_tests.cpp"            # ACTIVE - histogram and frequency analysis
    "algorithms/surface_integration_tests.cpp"  # ACTIVE - migrated to precision-aware
    "algorithms/volume_integration_tests.cpp"   # ACTIVE - migrated to precision-aware
    "algorithms/fourier_dft_tests.cpp"          # ACTIVE - Fourier DFT implementation
    "algorithms/fourier_fft_tests.cpp"          # ACTIVE - Fourier FFT implementation
    "algorithms/fourier_realfft_tests.cpp"      # ACTIVE - Fourier Real FFT implementation
    "algorithms/fourier_convolution_tests.cpp"  # ACTIVE - Fourier Convolution via FFT
    "algorithms/fourier_correlation_tests.cpp"  # ACTIVE - Fourier Correlation via FFT
    "algorithms/fourier_windowing_tests.cpp"    # ACTIVE - Fourier Window Functions
    "algorithms/fourier_spectrum_tests.cpp"     # ACTIVE - Power Spectrum & Spectrogram
    "algorithms/fourier_series_tests.cpp"       # ACTIVE - Fourier Series decomposition
    "algorithms/fourier_basis_tests.cpp"        # ACTIVE - Fourier orthogonal basis
    "algorithms/fourier_tests.cpp"              # ACTIVE - Comprehensive Fourier test suite
    
    "algorithms/eigensolver_building_blocks_tests.cpp" # ACTIVE - migrated to precision-aware
    "algorithms/optimization_framework_tests.cpp"      # ACTIVE - Optimization framework tests
    "algorithms/bvp_shooting_method_tests.cpp"         # ACTIVE - BVP shooting method solver
    "algorithms/Statistics/datadescriptors_tests.cpp"  # ACTIVE - typed data statistics
)

file(GLOB PRECISION_TESTS
    # "precision/derivation_precision_tests.cpp"    # DISABLED - testing float with vector only
    # "precision/integration_precision_tests.cpp"   # DISABLED - testing float with vector only
)

file(GLOB TOOLS_TESTS
    "tools/dataloader_tests.cpp"        # ACTIVE - data loading infrastructure
)

file(GLOB TESTBED_TESTS
    "testbeds/statisticstestbed_tests.cpp"  # ACTIVE - reference datasets
)

# Infrastructure tests
set(INFRASTRUCTURE_TESTS
    "test_precision_helpers_tests.cpp"  # ACTIVE - test infrastructure verification
    "test_matchers_tests.cpp"           # ACTIVE - custom Catch2 matchers for Real types
    # "test_simple_matcher.cpp"           # DEBUG - simple test
)

add_executable(MML_Tests 
    ${BASE_TESTS}
    ${CORE_TESTS}
    ${ALGORITHMS_TESTS}
    ${PRECISION_TESTS}
    ${TOOLS_TESTS}
    ${TESTBED_TESTS}
    ${INFRASTRUCTURE_TESTS}
)

# Link against Catch2WithMain - provides default main()
target_link_libraries(MML_Tests PRIVATE Catch2::Catch2WithMain)
target_include_directories(MML_Tests PRIVATE ${CMAKE_SOURCE_DIR}/mml)

# Add /FS flag for MSVC to allow parallel compilation with shared PDB
if(MSVC)
    target_compile_options(MML_Tests PRIVATE /FS)
endif()

# Use catch_discover_tests for better CTest integration
include(Catch)
catch_discover_tests(MML_Tests)

# ============================================================================
# Single Header Test - DISABLED pending generator fixes
# ============================================================================
# Issues found in MML.h generator (MinimalMathLibrary-m8rl):
# 1. File ordering: MMLPrecision.h must come before MMLBase.h
# 2. BOM characters in some source files not stripped
# 3. Missing <filesystem> header
# 4. Class redefinitions (ScalarFieldAnalyzer, VectorFieldAnalyzer)
# 5. Forward reference issues (DCT, EigenSolver, etc.)
#
# Uncomment when generator is fixed:
# add_executable(MML_SingleHeader_Tests single_header/test_single_header.cpp)
# target_link_libraries(MML_SingleHeader_Tests PRIVATE Catch2::Catch2WithMain)
# target_include_directories(MML_SingleHeader_Tests PRIVATE ${CMAKE_SOURCE_DIR}/mml)
# catch_discover_tests(MML_SingleHeader_Tests)
