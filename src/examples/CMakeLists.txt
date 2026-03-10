project(MML_Examples VERSION 1.0)

# Example 00: N-body gravity simulation (PRIME SHOWCASE)
add_executable(Example00_NBodyGravity
    00_N_body_gravity/main.cpp
)
target_include_directories(Example00_NBodyGravity PRIVATE 
    ${CMAKE_SOURCE_DIR}/mml
)

# Example 01: Projectile launch
add_executable(Example01_ProjectileLaunch
    01_projectile_launch/main.cpp
)
target_include_directories(Example01_ProjectileLaunch PRIVATE 
    ${CMAKE_SOURCE_DIR}/mml
)

# Example 02: Double pendulum
add_executable(Example02_DoublePendulum
    02_double_pendulum/main.cpp
)
target_include_directories(Example02_DoublePendulum PRIVATE 
    ${CMAKE_SOURCE_DIR}/mml
)

# Example 03: Formula 1 G-Force analysis
add_executable(Example03_F1GForces
    03_formula_1_sim/main.cpp
)
target_include_directories(Example03_F1GForces PRIVATE ${CMAKE_SOURCE_DIR}/mml)

# Example 04: 2D collision simulator
add_executable(Example04_CollisionSim2D
    04_collision_simulator_2d/main.cpp
)
target_include_directories(Example04_CollisionSim2D PRIVATE ${CMAKE_SOURCE_DIR}/mml)

# Example 05: Rigid body dynamics (self-contained - no MPL dependency)
add_executable(Example05_RigidBody
    05_rigid_body/main.cpp
)
target_include_directories(Example05_RigidBody PRIVATE 
    ${CMAKE_SOURCE_DIR}/mml
    ${CMAKE_CURRENT_SOURCE_DIR}/05_rigid_body
)

# Example 06: Lorentz transformations
add_executable(Example06_LorentzTransform
    06_Lorentz_transformations/main.cpp
)
target_include_directories(Example06_LorentzTransform PRIVATE ${CMAKE_SOURCE_DIR}/mml)

# Example 07: EM fields (TODO)
# add_executable(Example07_EMFields
#     07_EM_fields/main.cpp
# )
# target_include_directories(Example07_EMFields PRIVATE ${CMAKE_SOURCE_DIR}/mml)
