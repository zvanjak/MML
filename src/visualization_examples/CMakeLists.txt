project(MML_VisualizationApp VERSION 1.0)

# Collect all visualization example source files (excluding standalone utilities)
file(GLOB VISUALIZATION_SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp"
)
# Remove standalone utility programs from the main app
list(FILTER VISUALIZATION_SOURCES EXCLUDE REGEX "generate_.*\\.cpp$")

add_executable(MML_VisualizationApp 
    ${VISUALIZATION_SOURCES}
)

target_include_directories(MML_VisualizationApp PRIVATE ${CMAKE_SOURCE_DIR}/mml)

# Standalone utility: Generate 3D scalar function data for visualization testing
add_executable(Generate_Scalar3D_Data
    ${CMAKE_CURRENT_SOURCE_DIR}/generate_scalar3d_data.cpp
)
target_include_directories(Generate_Scalar3D_Data PRIVATE ${CMAKE_SOURCE_DIR}/mml)

# Standalone utility: Generate parametric surface data for visualization testing
add_executable(Generate_ParametricSurface_Data
    ${CMAKE_CURRENT_SOURCE_DIR}/generate_parametric_surface_data.cpp
)
target_include_directories(Generate_ParametricSurface_Data PRIVATE ${CMAKE_SOURCE_DIR}/mml)

# Copy visualization data files to build directory for runtime access
# Use GLOB_RECURSE to get files, and copy preserving directory structure
file(GLOB_RECURSE VISUALIZATION_DATA_FILES 
    "${CMAKE_CURRENT_SOURCE_DIR}/visualization_data/*.mml"
    "${CMAKE_CURRENT_SOURCE_DIR}/visualization_data/*.json"
    "${CMAKE_CURRENT_SOURCE_DIR}/visualization_data/*.csv"
    "${CMAKE_CURRENT_SOURCE_DIR}/visualization_data/*.txt"
)
foreach(DATA_FILE ${VISUALIZATION_DATA_FILES})
    file(RELATIVE_PATH REL_PATH "${CMAKE_CURRENT_SOURCE_DIR}/visualization_data" ${DATA_FILE})
    get_filename_component(REL_DIR ${REL_PATH} DIRECTORY)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/visualization_data/${REL_DIR}")
    configure_file(${DATA_FILE} "${CMAKE_CURRENT_BINARY_DIR}/visualization_data/${REL_PATH}" COPYONLY)
endforeach()
