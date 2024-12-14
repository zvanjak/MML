cmake_minimum_required(VERSION 3.20.0)

set (CMAKE_CXX_STANDARD 20)

project(MML_Demo VERSION 1.0)

add_compile_definitions(MML_USE_SINGLE_HEADER)

project(MML_DemoApp VERSION 1.0)
set(MML_DEMO_APP_SOURCES  src/demo_app/demo_app_main.cpp 

                          src/demo_app/base/demo_baseutils.cpp
                          src/demo_app/base/demo_geometry.cpp
                          src/demo_app/base/demo_matrix.cpp
                          src/demo_app/base/demo_matrix_other.cpp
                          src/demo_app/base/demo_matrixnm.cpp
                          src/demo_app/base/demo_tensors.cpp
                          src/demo_app/base/demo_vector.cpp
                          src/demo_app/base/demo_vectorn.cpp
                          src/demo_app/base/demo_geometry_2d.cpp
                          src/demo_app/base/demo_geometry_3d.cpp
                          src/demo_app/base/demo_polynom.cpp

                          src/demo_app/core/demo_coord_transf.cpp
                          src/demo_app/core/demo_coord_system.cpp
                          src/demo_app/core/demo_curves_surfaces.cpp
                          src/demo_app/core/demo_derivation.cpp
                          src/demo_app/core/demo_function.cpp
                          src/demo_app/core/demo_integration.cpp
                          src/demo_app/core/demo_interpolated_function.cpp
                          src/demo_app/core/demo_linear_alg_eq_solvers.cpp
                          src/demo_app/core/demo_serialization.cpp
                          src/demo_app/core/demo_metric_tensor.cpp
                          src/demo_app/core/demo_field_operations.cpp
                          src/demo_app/core/demo_coord_system_old.cpp
                          src/demo_app/core/demo_coord_transf_vector.cpp                                
                          src/demo_app/core/demo_dirac_function.cpp
                          src/demo_app/core/demo_fields.cpp
                          
                          src/demo_app/core/demo_intervals.cpp

                          src/demo_app/algorithms/demo_diff_geometry.cpp
                                                          
                          src/demo_app/algorithms/demo_function_analyzer.cpp
                          src/demo_app/algorithms/demo_ode_system_solvers.cpp
                          src/demo_app/algorithms/demo_path_integration.cpp
                          src/demo_app/algorithms/demo_root_finding.cpp
                          src/demo_app/algorithms/demo_statistics.cpp
                          src/demo_app/algorithms/demo_surface_integration.cpp
                          src/demo_app/algorithms/demo_volume_integration.cpp

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
                      )
add_executable(MML_DemoApp ${MML_DEMO_APP_SOURCES})
target_include_directories(MML_DemoApp PRIVATE include)


project(DocsApp VERSION 1.0)
set(MML_DOCS_APP_SOURCES    src/docs_demos/docs_app_main.cpp

                            src/docs_demos/docs_demo_coord_transf.cpp
                            src/docs_demos/docs_demo_derivation.cpp
                            src/docs_demos/docs_demo_diff_geometry.cpp
                            src/docs_demos/docs_demo_field_operations.cpp
                            src/docs_demos/docs_demo_functions.cpp
                            src/docs_demos/docs_demo_integration.cpp
                            src/docs_demos/docs_demo_integration_multidim.cpp
                            src/docs_demos/docs_demo_integration_path.cpp
                            src/docs_demos/docs_demo_integration_surface.cpp
                            src/docs_demos/docs_demo_integration_volume.cpp
                            src/docs_demos/docs_demo_interpolated_functions.cpp
                            src/docs_demos/docs_demo_intervals.cpp
                            src/docs_demos/docs_demo_geometry_2d_3d.cpp
                            src/docs_demos/docs_demo_lin_alg_solvers.cpp
                            src/docs_demos/docs_demo_matrix.cpp
                            src/docs_demos/docs_demo_matrixnm.cpp
                            src/docs_demos/docs_demo_ode_solvers.cpp
                            src/docs_demos/docs_demo_polynom.cpp
                            src/docs_demos/docs_demo_root_finding.cpp
                            src/docs_demos/docs_demo_tensors.cpp
                            src/docs_demos/docs_demo_vector.cpp
                            src/docs_demos/docs_demo_vectorn.cpp
                            src/docs_demos/docs_demo_visualizers.cpp
)
add_executable(DocsApp ${MML_DOCS_APP_SOURCES} )
target_include_directories(DocsApp PRIVATE include)


project(ReadmeApp VERSION 1.0)
set(MML_README_APP_SOURCES  src/readme_examples/readme_examples_main.cpp

                            src/readme_examples/readme1_vectors.cpp
                            src/readme_examples/readme2_matrices.cpp
                            src/readme_examples/readme3_tensors.cpp
                            src/readme_examples/readme4_geometry.cpp
                            src/readme_examples/readme5_polynoms.cpp
                            src/readme_examples/readme6_functions.cpp
                            src/readme_examples/readme7_derivation.cpp
                            src/readme_examples/readme8_integration.cpp
                            src/readme_examples/readme9_vector_field_operations.cpp
                            src/readme_examples/readme10_curves_surfaces.cpp
                            src/readme_examples/readme11_coord_system_transf.cpp
                            src/readme_examples/readme12_lin_alg_solvers.cpp
                            src/readme_examples/readme13_ode_solvers.cpp
                            src/readme_examples/readme14_function_analyzer.cpp
                            src/readme_examples/readme15_visualizators.cpp
)
add_executable(ReadmeApp ${MML_README_APP_SOURCES} )
target_include_directories(ReadmeApp PRIVATE include)


project(ExamplesApp VERSION 1.0)
set(MML_EXAMPLES_APP_SOURCES    src/examples/examples_main.cpp

                            src/examples/example_main.cpp
                            src/examples/example1_kosi_hitac.cpp
                            src/examples/example2_collision_calculator.cpp
                            src/examples/example3_tensor_of_inertia.cpp
                            src/examples/example4_gravity_field_investigations.cpp
                            src/examples/example5_voyager_travels.cpp
                            src/examples/example6_electric_charge_distribution.cpp
                            src/examples/example7_covariant_derivation.cpp
                            src/examples/example8_EM_field_investigations.cpp
                            src/examples/example9_thermodynamics_simulator.cpp
)
add_executable(ExamplesApp ${MML_EXAMPLES_APP_SOURCES} )
target_include_directories(ExamplesApp PRIVATE include)

 
set(CATCH_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/catch)
add_library(Catch INTERFACE)
target_include_directories(Catch INTERFACE ${CATCH_INCLUDE_DIR})

file(GLOB MML_TESTS tests/base/baseutils_tests.cpp
                    tests/base/geometry_tests.cpp
                    tests/base/geometry_2d_tests.cpp
                    tests/base/geometry_3d_tests.cpp
                    tests/base/intervals_tests.cpp
                    tests/base/matrix_tests.cpp
                    tests/base/matrix_sym_tests.cpp
                    tests/base/matrix_tridiag_tests.cpp
                    tests/base/matrixNM_tests.cpp
                    tests/base/tensors_tests.cpp
                    tests/base/vector_tests.cpp
                    tests/base/vectorN_tests.cpp 
                    tests/base/vectorN_subtypes_tests.cpp 
                    tests/base/polynom_tests.cpp

                    tests/core/coord_system_tests.cpp
                    tests/core/coord_transf_tests.cpp
                    tests/core/coord_transf_vector_tests.cpp
                    tests/core/coord_transf_tensor_tests.cpp
                    tests/core/derivation_tests.cpp                    
                    tests/core/field_operations_tests.cpp 
                    tests/core/function_tests.cpp
                    tests/core/integration_tests.cpp
                    tests/core/interpolated_functions_tests.cpp
                    tests/core/linear_alg_eq_solvers_tests.cpp 
                    tests/core/matrixutils_tests.cpp
                    tests/core/metric_tensor_tests.cpp

                    tests/algorithms/function_analyzer_tests.cpp
                    tests/algorithms/ode_system_solvers_tests.cpp  
                    tests/algorithms/parametric_curve_analyzer_tests.cpp  
                    tests/algorithms/path_integration_tests.cpp  
                    tests/algorithms/root_finding_tests.cpp  
                    tests/algorithms/surface_integration_tests.cpp  
                    tests/algorithms/volume_integration_tests.cpp  

                    tests/precision/derivation_precision_tests.cpp 
                    tests/precision/integration_precision_tests.cpp 
                )

add_executable(MML_Tests "tests/mml_test_runner.cpp" ${MML_TESTS} )
target_link_libraries(MML_Tests Catch)
target_include_directories(MML_Tests PRIVATE include)
