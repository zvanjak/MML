project(MML_VisualizationApp VERSION 1.0)

# Collect all visualization example source files
file(GLOB VISUALIZATION_SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp"
)

add_executable(MML_VisualizationApp 
    ${VISUALIZATION_SOURCES}
)

target_include_directories(MML_VisualizationApp PRIVATE ${CMAKE_SOURCE_DIR}/mml)

# Copy visualization data files to build directory for runtime access
file(GLOB VISUALIZATION_DATA_FILES "${CMAKE_CURRENT_SOURCE_DIR}/visualization_data/*")
foreach(DATA_FILE ${VISUALIZATION_DATA_FILES})
    get_filename_component(FILE_NAME ${DATA_FILE} NAME)
    configure_file(${DATA_FILE} ${CMAKE_CURRENT_BINARY_DIR}/visualization_data/${FILE_NAME} COPYONLY)
endforeach()
