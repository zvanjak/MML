project(MML_Examples VERSION 1.0)

# Example 00: N-body gravity simulation (PRIME SHOWCASE)
add_executable(Example00_NBodyGravity
    00_N_body_gravity/main.cpp
)
target_include_directories(Example00_NBodyGravity PRIVATE 
    ${CMAKE_SOURCE_DIR}/mml
    ${CMAKE_SOURCE_DIR}/src/book
)

# Example 01: Projectile launch
add_executable(Example01_ProjectileLaunch
    01_projectile_launch/main.cpp
)
target_include_directories(Example01_ProjectileLaunch PRIVATE 
    ${CMAKE_SOURCE_DIR}/mml
    ${CMAKE_SOURCE_DIR}/src/book
)

# Example 02: Double pendulum
add_executable(Example02_DoublePendulum
    02_double_pendulum/main.cpp
)
target_include_directories(Example02_DoublePendulum PRIVATE 
    ${CMAKE_SOURCE_DIR}/mml
    ${CMAKE_SOURCE_DIR}/src/book
)

# Example 03: 2D collision simulator
add_executable(Example03_CollisionSim2D
    03_collision_simulator_2d/main.cpp
)
target_include_directories(Example03_CollisionSim2D PRIVATE ${CMAKE_SOURCE_DIR}/mml)

# Example 04: Lorentz transformations
add_executable(Example04_LorentzTransform
    04_Lorentz_transformations/main.cpp
)
target_include_directories(Example04_LorentzTransform PRIVATE ${CMAKE_SOURCE_DIR}/mml)
