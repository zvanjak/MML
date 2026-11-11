project(MML_Examples VERSION 1.0)

# Example 00: N-body gravity simulation (PRIME SHOWCASE)
add_executable(Example00_NBodyGravity
    00_N_body_gravity/main.cpp
)
target_include_directories(Example00_NBodyGravity PRIVATE 
    ${CMAKE_SOURCE_DIR}
)

# Example 01: Projectile launch
add_executable(Example01_ProjectileLaunch
    01_projectile_launch/main.cpp
)
target_include_directories(Example01_ProjectileLaunch PRIVATE 
    ${CMAKE_SOURCE_DIR}
)

# Example 02: Double pendulum
add_executable(Example02_DoublePendulum
    02_double_pendulum/main.cpp
)
target_include_directories(Example02_DoublePendulum PRIVATE 
    ${CMAKE_SOURCE_DIR}
)

# Example 03: Formula 1 G-Force analysis
add_executable(Example03_F1GForces
    03_formula_1_sim/main.cpp
)
target_include_directories(Example03_F1GForces PRIVATE ${CMAKE_SOURCE_DIR})

# Example 04: 2D collision simulator
add_executable(Example04_CollisionSim2D
    04_collision_simulator_2d/main.cpp
)
target_include_directories(Example04_CollisionSim2D PRIVATE ${CMAKE_SOURCE_DIR})

# Example 05: Rigid body dynamics (self-contained - no MPL dependency)
add_executable(Example05_RigidBody
    05_rigid_body/main.cpp
)
target_include_directories(Example05_RigidBody PRIVATE 
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/05_rigid_body
)

# Example 06: Lorentz transformations (self-contained - no MPL dependency)
add_executable(Example06_LorentzTransform
    06_Lorentz_transformations/main.cpp
)
target_include_directories(Example06_LorentzTransform PRIVATE ${CMAKE_SOURCE_DIR})

# Example 07 requires MPL headers not distributed with public MML.
message(STATUS "Skipping Example07_EMFields: MPL dependencies are not part of public MML")

# Example 08: General relativity curvature and geodesics
add_executable(Example08_GRExamples
    08_GR_examples/main.cpp
)
target_include_directories(Example08_GRExamples PRIVATE ${CMAKE_SOURCE_DIR})

# Example 09: Typed differential forms gradient and directional derivative
add_executable(Example09_TypedFormsGradient
    09_typed_forms_gradient/main.cpp
)
target_include_directories(Example09_TypedFormsGradient PRIVATE ${CMAKE_SOURCE_DIR})

# Example 10: Typed differential forms flux and surface normal
add_executable(Example10_TypedFormsFlux
    10_typed_forms_flux/main.cpp
)
target_include_directories(Example10_TypedFormsFlux PRIVATE ${CMAKE_SOURCE_DIR})

# Example 11: Vector-space semantic layer
add_executable(Example11_VectorSpaces
    11_vector_spaces/main.cpp
)
target_include_directories(Example11_VectorSpaces PRIVATE ${CMAKE_SOURCE_DIR})

# Example 12: Tourism dataset statistics
add_executable(Example12_Statistics
    12_statistics/main.cpp
    12_statistics/TourismData.cpp
)
target_link_libraries(Example12_Statistics PRIVATE minimalmathlib::minimalmathlib)
if(MSVC)
    target_compile_options(Example12_Statistics PRIVATE /utf-8)
endif()
