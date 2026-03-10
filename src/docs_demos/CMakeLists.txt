project(MML_DocsApp VERSION 1.0)

# Collect sources from subdirectories (excluding mml_packages-dependent: optimization/, pde/)
file(GLOB BASE_SOURCES "base/*.cpp")
file(GLOB CORE_SOURCES "core/*.cpp")
file(GLOB ALGO_SOURCES "algorithms/*.cpp")
file(GLOB TOOLS_SOURCES "tools/*.cpp")
file(GLOB SYSTEMS_SOURCES "systems/*.cpp")
file(GLOB MAIN_SOURCES "*.cpp")

# Exclude verification files from MML_DocsApp (they have their own main() functions)
list(FILTER MAIN_SOURCES EXCLUDE REGEX "api_cheatsheet\\.cpp$")
list(FILTER MAIN_SOURCES EXCLUDE REGEX "quick_start_guide\\.cpp$")

# Exclude geometry demo - uses non-existent APIs (Point2Cartesian ostream, CoordTransf functions, etc)
list(FILTER BASE_SOURCES EXCLUDE REGEX "docs_demo_geometry\\.cpp$")

# Exclude demos that depend on mml_packages (fourier, optimization)
list(FILTER ALGO_SOURCES EXCLUDE REGEX "docs_demo_fourier\\.cpp$")
list(FILTER ALGO_SOURCES EXCLUDE REGEX "docs_demo_optimization\\.cpp$")

add_executable(MML_DocsApp 
    ${MAIN_SOURCES}
    ${BASE_SOURCES}
    ${CORE_SOURCES}
    ${ALGO_SOURCES}
    ${TOOLS_SOURCES}
    ${SYSTEMS_SOURCES}
)

target_include_directories(MML_DocsApp PRIVATE 
    ${CMAKE_SOURCE_DIR}/mml
    ${CMAKE_SOURCE_DIR}/test_data
    ${CMAKE_SOURCE_DIR}/test_beds
)

# Organize in IDE (Visual Studio folders)
source_group("Main" FILES ${MAIN_SOURCES})
source_group("Base" FILES ${BASE_SOURCES})
source_group("Core" FILES ${CORE_SOURCES})
source_group("Algorithms" FILES ${ALGO_SOURCES})
source_group("Tools" FILES ${TOOLS_SOURCES})
source_group("Systems" FILES ${SYSTEMS_SOURCES})

# Documentation guide verification executables
add_executable(API_Cheatsheet_Verification api_cheatsheet.cpp)
target_include_directories(API_Cheatsheet_Verification PRIVATE ${CMAKE_SOURCE_DIR}/mml)

add_executable(Quick_Start_Guide_Verification quick_start_guide.cpp)
target_include_directories(Quick_Start_Guide_Verification PRIVATE ${CMAKE_SOURCE_DIR}/mml)
