# CMakeLists.txt for MPL (Minimal Physics Library) Tests
# 
# Unit tests for the physics simulation components in MPL.
# Uses Catch2 for testing framework (same as main MML tests).
#
# Test Categories:
# - SpecialRelativity: Lorentz transformations, relativistic mechanics
# - RigidBody: Rigid body dynamics, moment of inertia, collisions
# - Electromagnetism: EM tensor, field calculations
# - Mechanics: Projectiles, pendulums, oscillators
# - Gravity: N-body, orbital mechanics

project(MPL_Tests VERSION 1.0)

# Catch2 is fetched in root CMakeLists.txt via FetchContent

###############################################################################
# SPECIAL RELATIVITY TESTS
###############################################################################
set(SPECIAL_RELATIVITY_TESTS
    "special_relativity/lorentz_transformation_tests.cpp"
    "special_relativity/time_dilation_tests.cpp"
    "special_relativity/length_contraction_tests.cpp"
)

###############################################################################
# RIGID BODY TESTS
###############################################################################
set(RIGID_BODY_TESTS
    "rigid_body/moment_of_inertia_tests.cpp"
    "rigid_body/rigid_body_state_tests.cpp"
    "rigid_body/euler_equations_tests.cpp"
    "rigid_body/sphere_tests.cpp"
    "rigid_body/box_inertia_tests.cpp"
)

###############################################################################
# ELECTROMAGNETISM TESTS
###############################################################################
set(ELECTROMAGNETISM_TESTS
    "electromagnetism/em_tensor_tests.cpp"
    "electromagnetism/field_transformation_tests.cpp"
)

###############################################################################
# MECHANICS TESTS (future)
###############################################################################
set(MECHANICS_TESTS
    # "mechanics/projectile_tests.cpp"
    # "mechanics/pendulum_tests.cpp"
)

###############################################################################
# ALL MPL TESTS
###############################################################################
set(ALL_MPL_TESTS
    ${SPECIAL_RELATIVITY_TESTS}
    ${RIGID_BODY_TESTS}
    ${ELECTROMAGNETISM_TESTS}
    ${MECHANICS_TESTS}
    mpl_test_main.cpp
)

# Create test executable
add_executable(MPL_Tests ${ALL_MPL_TESTS})

target_include_directories(MPL_Tests PRIVATE
    ${CMAKE_SOURCE_DIR}/mml
    ${CMAKE_SOURCE_DIR}/src/book/mpl
    ${CMAKE_SOURCE_DIR}/tests          # For TestPrecision.h and TestMatchers.h
)

target_link_libraries(MPL_Tests PRIVATE Catch2::Catch2WithMain)

# Register with CTest
include(CTest)
include(Catch)
catch_discover_tests(MPL_Tests
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

# Add test label for filtering
set_tests_properties(${TESTS} PROPERTIES LABELS "mpl")
