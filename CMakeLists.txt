cmake_minimum_required(VERSION 3.0.0)

set (CMAKE_CXX_STANDARD 20)

project(MML_Demo VERSION 0.1.0)

set( MML_DEMO_SOURCES   src/main.cpp 

                        src/utilities/demo_intervals.cpp

                        src/core/demo_algebra.cpp
                        src/core/demo_coord_transf.cpp
                        src/core/demo_coreutils.cpp
                        src/core/demo_derivation.cpp
                        src/core/demo_function.cpp
                        src/core/demo_geometry.cpp
                        src/core/demo_integration.cpp
                        src/core/demo_interpolated_function.cpp
                        src/core/demo_linear_alg_eq_solvers.cpp
                        src/core/demo_linear_functional.cpp
                        src/core/demo_linear_operator.cpp
                        src/core/demo_matrix.cpp
                        src/core/demo_matrix_other.cpp
                        src/core/demo_matrixnm.cpp
                        src/core/demo_polynom.cpp
                        src/core/demo_quadratic_form.cpp
                        src/core/demo_tensors.cpp
                        src/core/demo_vector.cpp
                        src/core/demo_vectorn.cpp
                        src/core/demo_vector_spaces.cpp
                        src/core/demo_metric_tensor.cpp
                        src/core/demo_field_operations.cpp

                        src/basic_types/demo_coord_system.cpp
                        src/basic_types/demo_covar_contravar_transf.cpp                                
                        src/basic_types/demo_dirac_function.cpp
                        src/basic_types/demo_fields.cpp
                        src/basic_types/demo_function_space.cpp
                        src/basic_types/demo_geometry_2d.cpp
                        src/basic_types/demo_geometry_3d.cpp
                        src/basic_types/demo_interpolators.cpp
                        src/basic_types/demo_surfaces.cpp

                        src/algorithms/demo_diff_geometry.cpp
                        src/algorithms/demo_eigen_solvers.cpp                                
                        src/algorithms/demo_function_analyzer.cpp
                        src/algorithms/demo_ode_system_solvers.cpp
                        src/algorithms/demo_path_integration.cpp
                        src/algorithms/demo_root_finding.cpp
                        src/algorithms/demo_statistics.cpp

                        src/testing_speed/test_speed_functions.cpp
                        src/testing_speed/test_speed_derivation.cpp
                        src/testing_speed/test_speed_linear_alg_eq_solvers.cpp

                        src/testing_precision/test_precision_derivation.cpp
                        src/testing_precision/test_precision_integration.cpp
                        src/testing_precision/test_precision_interpolation.cpp
                            
                        src/readme_examples.cpp
                    )

add_executable(MML_Demo ${MML_DEMO_SOURCES})

add_compile_definitions(MML_USE_SINGLE_HEADER)
target_include_directories(MML_Demo PRIVATE include)


set(CATCH_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/catch)
add_library(Catch INTERFACE)
target_include_directories(Catch INTERFACE ${CATCH_INCLUDE_DIR})

file(GLOB MML_TESTS tests/core/algebra_tests.cpp
                    tests/core/coord_transf_tests.cpp
                    tests/core/coreutils_tests.cpp
                    tests/core/derivation_tests.cpp                    
                    tests/core/function_tests.cpp
                    tests/core/functional_tests.cpp
                    tests/core/geometry_tests.cpp
                    tests/core/integration_tests.cpp
                    tests/core/interpolated_functions_tests.cpp
                    tests/core/linear_alg_eq_solvers_tests.cpp 
                    tests/core/matrix_tests.cpp
                    tests/core/matrix_other_tests.cpp
                    tests/core/matrixNM_tests.cpp
                    tests/core/polynom_tests.cpp
                    tests/core/tensors_tests.cpp
                    tests/core/vector_tests.cpp
                    tests/core/vectorN_tests.cpp 
                    tests/core/field_operations_tests.cpp 
                    tests/core/metric_tensor_tests.cpp

                    tests/basic_types/coord_system_tests.cpp
                    tests/basic_types/function_space_tests.cpp
                    tests/core/geometry_2d_tests.cpp
                    tests/core/geometry_3d_tests.cpp
                    tests/basic_types/surfaces_tests.cpp

                    tests/algorithms/diff_geom_alg_tests.cpp  
                    tests/algorithms/eigen_solvers_tests.cpp 
                    tests/algorithms/function_analyzer_tests.cpp
                    tests/algorithms/ode_system_solvers_tests.cpp  
                    tests/algorithms/path_integration_tests.cpp  

                    tests/precision/derivation_precision_tests.cpp 
                    tests/precision/integration_precision_tests.cpp 

                    tests/intervals_tests.cpp
                    )

add_executable(MML_Tests "mmlib_test_runner.cpp" ${MML_TESTS} )
target_link_libraries(MML_Tests Catch)
target_include_directories(MML_Tests PRIVATE include)