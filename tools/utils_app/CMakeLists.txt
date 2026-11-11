# MML Utilities Application
# Part of the main MML build - do not build standalone

add_executable(MML_UtilsApp utils_app.cpp)

target_include_directories(MML_UtilsApp PRIVATE 
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}
)
