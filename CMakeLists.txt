cmake_minimum_required(VERSION 3.20.0)

set (CMAKE_CXX_STANDARD 20)

project(MML_Demo VERSION 1.0)

add_compile_definitions(MML_USE_SINGLE_HEADER)

set( MML_DEMO_SOURCES   src/demo_app_main.cpp 

                        src/base/demo_algebra.cpp
                        src/base/demo_baseutils.cpp
                        src/base/demo_geometry.cpp
                        src/base/demo_matrix.cpp
                        src/base/demo_matrix_other.cpp
                        src/base/demo_matrixnm.cpp
                        src/base/demo_tensors.cpp
                        src/base/demo_vector.cpp
                        src/base/demo_vectorn.cpp
                        src/base/demo_geometry_2d.cpp
                        src/base/demo_geometry_3d.cpp
                        src/base/demo_linear_functional.cpp
                        src/base/demo_linear_operator.cpp
                        src/base/demo_polynom.cpp
                        src/base/demo_quadratic_form.cpp
                        src/base/demo_vector_spaces.cpp

                        src/core/demo_chebyshev_approximation.cpp
                        src/core/demo_coord_transf.cpp
                        src/core/demo_coord_system.cpp
                        src/core/demo_derivation.cpp
                        src/core/demo_function.cpp
                        src/core/demo_integration.cpp
                        src/core/demo_interpolated_function.cpp
                        src/core/demo_linear_alg_eq_solvers.cpp
                        src/core/demo_serialization.cpp
                        src/core/demo_metric_tensor.cpp
                        src/core/demo_field_operations.cpp
                        src/core/demo_coord_system_old.cpp
                        src/core/demo_covar_contravar_transf.cpp                                
                        src/core/demo_dirac_function.cpp
                        src/core/demo_fields.cpp
                        src/core/demo_function_space.cpp
                        src/core/demo_interpolators.cpp
                        src/core/demo_surfaces.cpp
                        src/core/demo_intervals.cpp

                        src/algorithms/demo_diff_geometry.cpp
                        src/algorithms/demo_eigen_solvers.cpp                                
                        src/algorithms/demo_function_analyzer.cpp
                        src/algorithms/demo_ode_system_solvers.cpp
                        src/algorithms/demo_path_integration.cpp
                        src/algorithms/demo_root_finding.cpp
                        src/algorithms/demo_statistics.cpp
                        src/algorithms/demo_surface_integration.cpp
                        src/algorithms/demo_volume_integration.cpp

                        src/testing_speed/test_speed_functions.cpp
                        src/testing_speed/test_speed_derivation.cpp
                        src/testing_speed/test_speed_linear_alg_eq_solvers.cpp

                        src/testing_precision/test_precision_derivation.cpp
                        src/testing_precision/test_precision_integration.cpp
                        src/testing_precision/test_precision_interpolation.cpp

                        src/visualization_examples/demo1_lorenz_multi_func.cpp
                        src/visualization_examples/demo2_lorenz_parametric_curve.cpp
                        src/visualization_examples/demo3_surface.cpp
                        src/visualization_examples/demo4_vector_field_gravity.cpp
                        src/visualization_examples/demo5_vector_field_EM.cpp
                        src/visualization_examples/visualization_examples.cpp
                            
                        docs_examples/readme_examples/readme1_vectors_matrices.cpp
                        docs_examples/readme_examples/readme2_lin_alg_solvers.cpp
                        docs_examples/readme_examples/readme3_defining_functions.cpp
                        docs_examples/readme_examples/readme4_working_with_functions.cpp
                        docs_examples/readme_examples/readme5_ode_solvers.cpp
                        docs_examples/readme_examples/readme6_vector_field_operations.cpp
                        docs_examples/readme_examples/readme7_parametric_curves.cpp
                        docs_examples/readme_examples/readme8_visualizators.cpp
                        docs_examples/readme_examples/readme9_function_analyzer.cpp
                        docs_examples/readme_examples/main_readme_examples.cpp
                      )


add_executable(MML_Demo ${MML_DEMO_SOURCES})
target_include_directories(MML_Demo PRIVATE include)


project(DocsApp VERSION 1.0)
set(MML_DOCS_APP_SOURCES    docs_examples/docs_app_main.cpp

                            docs_examples/readme_examples/readme1_vectors_matrices.cpp
                            docs_examples/readme_examples/readme2_lin_alg_solvers.cpp
                            docs_examples/readme_examples/readme3_defining_functions.cpp
                            docs_examples/readme_examples/readme4_working_with_functions.cpp
                            docs_examples/readme_examples/readme5_ode_solvers.cpp
                            docs_examples/readme_examples/readme6_vector_field_operations.cpp
                            docs_examples/readme_examples/readme7_parametric_curves.cpp
                            docs_examples/readme_examples/readme8_visualizators.cpp
                            docs_examples/readme_examples/readme9_function_analyzer.cpp
                            
                            docs_examples/readme_examples/main_readme_examples.cpp

                            docs_examples/examples/example_main.cpp
                            docs_examples/examples/example1_kosi_hitac.cpp
                            docs_examples/examples/example2_collision_calculator.cpp
                            docs_examples/examples/example3_tensor_of_inertia.cpp
                            docs_examples/examples/example4_gravity_field_investigations.cpp
                            docs_examples/examples/example5_voyager_travels.cpp
                            docs_examples/examples/example6_electric_charge_distribution.cpp
                            docs_examples/examples/example7_covariant_derivation.cpp
                            docs_examples/examples/example8_EM_field_investigations.cpp
                            docs_examples/examples/example9_thermodynamics_simulator.cpp
                            
                            docs_examples/examples/main_docs_examples.cpp

                            docs_examples/demos/docs_demo_coord_transf.cpp
                            docs_examples/demos/docs_demo_derivation.cpp
                            docs_examples/demos/docs_demo_diff_geometry.cpp
                            docs_examples/demos/docs_demo_eigen_solvers.cpp
                            docs_examples/demos/docs_demo_field_operations.cpp
                            docs_examples/demos/docs_demo_functions.cpp
                            docs_examples/demos/docs_demo_integration.cpp
                            docs_examples/demos/docs_demo_integration_multidim.cpp
                            docs_examples/demos/docs_demo_integration_path.cpp
                            docs_examples/demos/docs_demo_integration_surface.cpp
                            docs_examples/demos/docs_demo_interpolated_functions.cpp
                            docs_examples/demos/docs_demo_geometry_2d_3d.cpp
                            docs_examples/demos/docs_demo_lin_alg_solvers.cpp
                            docs_examples/demos/docs_demo_matrix.cpp
                            docs_examples/demos/docs_demo_matrixnm.cpp
                            docs_examples/demos/docs_demo_ode_solvers.cpp
                            docs_examples/demos/docs_demo_polynom.cpp
                            docs_examples/demos/docs_demo_root_finding.cpp
                            docs_examples/demos/docs_demo_surface_integration.cpp
                            docs_examples/demos/docs_demo_tensors.cpp
                            docs_examples/demos/docs_demo_vector.cpp
                            docs_examples/demos/docs_demo_visualizers.cpp
                            docs_examples/demos/docs_demo_volume_integration.cpp
                            
                            docs_examples/demos/main_docs_demos.cpp
)
add_executable(DocsApp ${MML_DOCS_APP_SOURCES})
target_include_directories(DocsApp PRIVATE include)


set(CATCH_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/catch)
add_library(Catch INTERFACE)
target_include_directories(Catch INTERFACE ${CATCH_INCLUDE_DIR})

file(GLOB MML_TESTS tests/base/algebra_tests.cpp

                    tests/base/baseutils_tests.cpp
                    tests/base/geometry_tests.cpp
                    tests/base/matrix_tests.cpp
                    tests/base/matrix_sym_tests.cpp
                    tests/base/matrix_tridiag_tests.cpp
                    tests/base/matrix_banddiag_tests.cpp
                    tests/base/matrixNM_tests.cpp
                    tests/base/tensors_tests.cpp
                    tests/base/vector_tests.cpp
                    tests/base/vectorN_tests.cpp 
                    tests/base/vector_types_tests.cpp 
                    tests/base/geometry_2d_tests.cpp
                    tests/base/geometry_3d_tests.cpp
                    tests/base/polynom_tests.cpp

                    tests/core/coord_transf_tests.cpp
                    tests/core/derivation_tests.cpp                    
                    tests/core/field_operations_tests.cpp 
                    tests/core/function_tests.cpp
                    tests/core/function_serialization_tests.cpp
                    tests/core/functional_tests.cpp
                    tests/core/integration_tests.cpp
                    tests/core/interpolated_functions_tests.cpp
                    tests/core/intervals_tests.cpp
                    tests/core/linear_alg_eq_solvers_tests.cpp 
                    tests/core/matrixutils_tests.cpp
                    tests/core/metric_tensor_tests.cpp
                    tests/core/coord_system_tests.cpp
                    tests/core/function_space_tests.cpp

                    tests/algorithms/diff_geom_alg_tests.cpp  
                    tests/algorithms/eigen_solvers_tests.cpp 
                    tests/algorithms/function_analyzer_tests.cpp
                    tests/algorithms/ode_system_solvers_tests.cpp  
                    tests/algorithms/path_integration_tests.cpp  
                    tests/algorithms/root_finding_tests.cpp  
                    tests/algorithms/root_finding_multidim_tests.cpp  

                    tests/precision/derivation_precision_tests.cpp 
                    tests/precision/integration_precision_tests.cpp 
)

add_executable(MML_Tests "tests/mml_test_runner.cpp" ${MML_TESTS} )
target_link_libraries(MML_Tests Catch)
target_include_directories(MML_Tests PRIVATE include)