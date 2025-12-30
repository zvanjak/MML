cmake_minimum_required(VERSION 3.10)

add_executable(test_linux_viz_standalone test_linux_viz_standalone.cpp)
target_include_directories(test_linux_viz_standalone PRIVATE 
    ${CMAKE_SOURCE_DIR}/mml
)
target_compile_features(test_linux_viz_standalone PRIVATE cxx_std_17)

add_executable(test_backend_switching test_backend_switching.cpp)
target_include_directories(test_backend_switching PRIVATE 
    ${CMAKE_SOURCE_DIR}/mml
)
target_compile_features(test_backend_switching PRIVATE cxx_std_17)

add_executable(test_all_visualizers test_all_visualizers.cpp)
target_include_directories(test_all_visualizers PRIVATE 
    ${CMAKE_SOURCE_DIR}/mml
)
target_compile_features(test_all_visualizers PRIVATE cxx_std_17)
