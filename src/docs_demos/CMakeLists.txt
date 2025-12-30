project(MML_DocsApp VERSION 1.0)

# Collect sources from all subdirectories
file(GLOB BASE_SOURCES "base/*.cpp")
file(GLOB CORE_SOURCES "core/*.cpp")
file(GLOB ALGO_SOURCES "algorithms/*.cpp")
file(GLOB TOOLS_SOURCES "tools/*.cpp")
file(GLOB MAIN_SOURCES "*.cpp")

add_executable(MML_DocsApp 
    ${MAIN_SOURCES}
    ${BASE_SOURCES}
    ${CORE_SOURCES}
    ${ALGO_SOURCES}
    ${TOOLS_SOURCES}
)

target_include_directories(MML_DocsApp PRIVATE ${CMAKE_SOURCE_DIR}/mml)

# Organize in IDE (Visual Studio folders)
source_group("Main" FILES ${MAIN_SOURCES})
source_group("Base" FILES ${BASE_SOURCES})
source_group("Core" FILES ${CORE_SOURCES})
source_group("Algorithms" FILES ${ALGO_SOURCES})
source_group("Tools" FILES ${TOOLS_SOURCES})
