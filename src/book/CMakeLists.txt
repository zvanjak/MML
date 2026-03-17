# CMakeLists.txt for MML_BookApp
# Unified book application that includes all chapter examples
#
# This builds a single executable that can run all chapter examples.
# Each chapter also has its own standalone executable via the chapters/ subdirectory.

project(MML_BookApp VERSION 1.0)

# Main book_main.cpp
set(BOOK_MAIN_SOURCE "book_main.cpp")

# Collect chapter sources (excluding main.cpp files from each chapter)
# We get all .cpp files then filter out main.cpp

# Helper function to collect sources excluding main.cpp files
# For most chapters, main_chapter_XX.cpp IS the content file (main() was separated to main.cpp)
# Exception: Chapter 17 has duplicate code in main_chapter_17.cpp that we need to exclude
function(collect_chapter_sources VAR_NAME CHAPTER_DIR)
    file(GLOB ALL_SRCS "${CHAPTER_DIR}/*.cpp")
    set(FILTERED_SRCS)
    foreach(SRC ${ALL_SRCS})
        get_filename_component(FNAME ${SRC} NAME)
        # Exclude main.cpp (standalone main functions)
        if(NOT FNAME STREQUAL "main.cpp")
            # Also exclude main_chapter_17.cpp specifically (contains duplicates of lorentz_demos.cpp)
            if(NOT FNAME STREQUAL "main_chapter_17.cpp")
                list(APPEND FILTERED_SRCS ${SRC})
            endif()
        endif()
    endforeach()
    set(${VAR_NAME} ${FILTERED_SRCS} PARENT_SCOPE)
endfunction()

# Collect sources from all chapters
collect_chapter_sources(CH01_SOURCES "chapters/Chapter_01_basic_objects")
collect_chapter_sources(CH02_SOURCES "chapters/Chapter_02_vizualization")
collect_chapter_sources(CH03_SOURCES "chapters/Chapter_03_basic_algorithms")
collect_chapter_sources(CH04_SOURCES "chapters/Chapter_04_collision_simulator")
collect_chapter_sources(CH05_SOURCES "chapters/Chapter_05_throwing_thing_in_the_air")
collect_chapter_sources(CH06_SOURCES "chapters/Chapter_06_pendulum")
collect_chapter_sources(CH07_SOURCES "chapters/Chapter_07_double_and_spherical_pendulum")
collect_chapter_sources(CH08_SOURCES "chapters/Chapter_08_gravity")
collect_chapter_sources(CH09_SOURCES "chapters/Chapter_09_gravity_N_body_problem")
collect_chapter_sources(CH11_SOURCES "chapters/Chapter_11_coordinate_transformations")
collect_chapter_sources(CH12_SOURCES "chapters/Chapter_12_carousel")
collect_chapter_sources(CH13_SOURCES "chapters/Chapter_13_projectile_launch")
collect_chapter_sources(CH14_SOURCES "chapters/Chapter_14_rotations_quaternions")
collect_chapter_sources(CH15_SOURCES "chapters/Chapter_15_rigid_body")
collect_chapter_sources(CH17_SOURCES "chapters/Chapter_17_lorentz_transformations")
collect_chapter_sources(CH18_SOURCES "chapters/Chapter_18_twin_paradox")
collect_chapter_sources(CH19_SOURCES "chapters/Chapter_19_static_EM_fields")
collect_chapter_sources(CH20_SOURCES "chapters/Chapter_20_dynamic_EM_fields")
collect_chapter_sources(CH21_SOURCES "chapters/Chapter_21_diff_geometry_curves_surfaces")
collect_chapter_sources(CH22_SOURCES "chapters/Chapter_22_general_relativity")

# Combine all sources
set(ALL_CHAPTER_SOURCES
    ${CH01_SOURCES}
    ${CH02_SOURCES}
    ${CH03_SOURCES}
    ${CH04_SOURCES}
    ${CH05_SOURCES}
    ${CH06_SOURCES}
    ${CH07_SOURCES}
    ${CH08_SOURCES}
    ${CH09_SOURCES}
    ${CH11_SOURCES}
    ${CH12_SOURCES}
    ${CH13_SOURCES}
    ${CH14_SOURCES}
    ${CH15_SOURCES}
    ${CH17_SOURCES}
    ${CH18_SOURCES}
    ${CH19_SOURCES}
    ${CH20_SOURCES}
    ${CH21_SOURCES}
    ${CH22_SOURCES}
)

# Create the unified book executable
add_executable(MML_BookApp 
    ${BOOK_MAIN_SOURCE}
    ${ALL_CHAPTER_SOURCES}
)

target_include_directories(MML_BookApp PRIVATE 
    ${CMAKE_SOURCE_DIR}/mml
    ${CMAKE_SOURCE_DIR}/test_data
    ${CMAKE_SOURCE_DIR}/test_beds
    ${CMAKE_CURRENT_SOURCE_DIR}/mpl
)

# Note: Individual chapter applications are built by src/book/chapters/
# which is already added at the root CMakeLists.txt level

# Add MPL unit tests
add_subdirectory(mpl_tests)
