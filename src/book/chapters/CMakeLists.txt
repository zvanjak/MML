# CMakeLists.txt for Book Chapters
# Each chapter is a standalone application for the book "Exploring Physics with C++"
# This master file builds all chapter applications

cmake_minimum_required(VERSION 3.20)

# Add each chapter subdirectory
# Chapters are built as independent executables

# Chapter 01: Basic Objects
add_subdirectory(Chapter_01_basic_objects)

# Chapter 02: Visualization
add_subdirectory(Chapter_02_vizualization)

# Chapter 03: Basic Algorithms
add_subdirectory(Chapter_03_basic_algorithms)

# Chapter 04: Collision Simulator
add_subdirectory(Chapter_04_collision_simulator)

# Chapter 05: Throwing Things in the Air
add_subdirectory(Chapter_05_throwing_thing_in_the_air)

# Chapter 06: Pendulum
add_subdirectory(Chapter_06_pendulum)

# Chapter 07: Double and Spherical Pendulum
add_subdirectory(Chapter_07_double_and_spherical_pendulum)

# Chapter 08: Gravity
add_subdirectory(Chapter_08_gravity)

# Chapter 09: N-Body Gravity Problem
add_subdirectory(Chapter_09_gravity_N_body_problem)

# Chapter 10: N-Body Applications (no CMakeLists yet)
# add_subdirectory(Chapter_10_N_body_applications)

# Chapter 11: Coordinate Transformations
add_subdirectory(Chapter_11_coordinate_transformations)

# Chapter 12: Carousel
add_subdirectory(Chapter_12_carousel)

# Chapter 13: Projectile Launch
add_subdirectory(Chapter_13_projectile_launch)

# Chapter 14: Rotations and Quaternions
add_subdirectory(Chapter_14_rotations_quaternions)

# Chapter 15: Rigid Body
add_subdirectory(Chapter_15_rigid_body)

# Chapter 16: Rigid Body Simulator (no CMakeLists yet)
# add_subdirectory(Chapter_16_rigid_body_simulator)

# Chapter 17: Lorentz Transformations
add_subdirectory(Chapter_17_lorentz_transformations)

# Chapter 18: Twin Paradox
add_subdirectory(Chapter_18_twin_paradox)

# Chapter 19: Static EM Fields
add_subdirectory(Chapter_19_static_EM_fields)

# Chapter 20: Dynamic EM Fields
add_subdirectory(Chapter_20_dynamic_EM_fields)

# Chapter 21: Differential Geometry (Curves and Surfaces)
add_subdirectory(Chapter_21_diff_geometry_curves_surfaces)

# Chapter 22: General Relativity
add_subdirectory(Chapter_22_general_relativity)
