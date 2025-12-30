project(MML_PrecisionTestingApp VERSION 1.0)

# Precision testing framework and test files
set(PRECISION_TEST_SOURCES
    main.cpp
    test_precision_derivation.cpp
    test_precision_integration.cpp
    test_precision_interpolation.cpp
    test_precision_ode.cpp
    test_precision_roots.cpp
    test_precision_linalg.cpp
)

# Note: test_framework_verify.cpp has its own main() - excluded from main app
# Can be built separately if needed

add_executable(MML_PrecisionTestingApp 
    ${PRECISION_TEST_SOURCES}
)

target_include_directories(MML_PrecisionTestingApp PRIVATE 
    ${CMAKE_SOURCE_DIR}/mml
    ${CMAKE_SOURCE_DIR}/test_data
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# Set C++17 standard
target_compile_features(MML_PrecisionTestingApp PRIVATE cxx_std_17)

# Standalone framework verification test (optional)
add_executable(MML_FrameworkVerify
    test_framework_verify.cpp
)
target_include_directories(MML_FrameworkVerify PRIVATE 
    ${CMAKE_SOURCE_DIR}/mml
    ${CMAKE_CURRENT_SOURCE_DIR}
)
target_compile_features(MML_FrameworkVerify PRIVATE cxx_std_17)
