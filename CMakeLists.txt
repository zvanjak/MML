cmake_minimum_required(VERSION 3.0.0)

set (CMAKE_CXX_STANDARD 20)

project(MML_Demo VERSION 0.1.0)

set( MML_DEMO_SOURCES   src/main.cpp 

                        src/demo_algebra.cpp
                        src/demo_coord_system.cpp
                        src/demo_coord_transf.cpp
                        src/demo_covar_contravar_transf.cpp                                
                        src/demo_derivation.cpp
                        src/demo_diff_eq_solvers.cpp
                        src/demo_eigen_solvers.cpp                                
                        src/demo_field_operations.cpp
                        src/demo_fields.cpp
                        src/demo_geometry.cpp
                        src/demo_integration.cpp
                        src/demo_interpolators.cpp
                        src/demo_linear_alg_eq_solvers.cpp
                        src/demo_path_integration.cpp
                        src/demo_diff_geometry.cpp

                        src/demo_type_function.cpp
                        src/demo_type_interpolated_function.cpp
                        src/demo_type_matrix.cpp
                        src/demo_type_matrixnm.cpp
                        src/demo_type_polynom.cpp
                        src/demo_type_tensors.cpp
                        src/demo_type_vector.cpp
                        src/demo_type_vectorn.cpp

                        src/testing_speed/test_speed_functions.cpp
                        src/testing_speed/test_speed_derivation.cpp
                        src/testing_speed/test_speed_linear_alg_eq_solvers.cpp

                        src/testing_precision/test_precision_derivation.cpp
                        src/testing_precision/test_precision_integration.cpp

                        src/readme_examples.cpp
                    )

add_executable(MML_Demo ${MML_DEMO_SOURCES})

add_compile_definitions(MML_USE_SINGLE_HEADER)
target_include_directories(MML_Demo PRIVATE include)

set(CATCH_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/catch)
add_library(Catch INTERFACE)
target_include_directories(Catch INTERFACE ${CATCH_INCLUDE_DIR})

file(GLOB MML_TESTS tests/derivation_tests.cpp 
                    tests/diff_eq_solvers_tests.cpp  
                    tests/diff_geom_alg_tests.cpp  
                    tests/eigen_solvers_tests.cpp 
                    tests/function_tests.cpp
                    tests/function_analyzer_tests.cpp
                    tests/geometry_tests.cpp
                    tests/integration_tests.cpp
                    tests/interpolated_functions_tests.cpp
                    tests/intervals_tests.cpp
                    tests/linear_alg_eq_solvers_tests.cpp 
                    tests/matrix_tests.cpp
                    tests/matrixNM_tests.cpp
                    tests/polynom_tests.cpp
                    tests/utilities_tests.cpp
                    tests/tensors_tests.cpp
                    tests/vector_tests.cpp
                    tests/vectorN_tests.cpp  
                    )

add_executable(MML_Tests "mmlib_test_runner.cpp" ${MML_TESTS} )
target_link_libraries(MML_Tests Catch)
target_include_directories(MML_Tests PRIVATE include)