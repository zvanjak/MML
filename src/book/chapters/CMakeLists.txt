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

# Chapter 10: Coordinate Transformations
add_subdirectory(Chapter_10_coordinate_transformations)

# Chapter 11: Carousel
add_subdirectory(Chapter_11_carousel)

# Chapter 12: Projectile Launch
add_subdirectory(Chapter_12_projectile_launch)

# Chapter 13: Rigid Body
add_subdirectory(Chapter_13_rigid_body)

# Chapter 14: Rotations and Quaternions
add_subdirectory(Chapter_14_rotations_quaternions)

# Chapter 15: Lorentz Transformations
add_subdirectory(Chapter_15_lorentz_transformations)

# Chapter 17: Static Electric Fields
add_subdirectory(Chapter_17_static_electric_fields)

# Chapter 18: Static Magnetic Fields
add_subdirectory(Chapter_18_static_magnetic_fields)

# Chapter 19: Dynamic EM Fields
add_subdirectory(Chapter_19_dynamic_EM_fields)

# Chapter 20: Differential Geometry (Curves and Surfaces)
add_subdirectory(Chapter_20_diff_geometry_curves_surfaces)

# Chapter 21: General Relativity
add_subdirectory(Chapter_21_general_relativity)

# Note: Chapter 16 (Twin Paradox) is optional and may stay integrated in Chapter 15
