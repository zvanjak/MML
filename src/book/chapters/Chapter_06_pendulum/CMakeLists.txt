project(Chapter06App VERSION 1.0)

set(CHAPTER_SOURCES
    main.cpp
    chapter_06_coordinator.cpp
    simple_pendulum.cpp
    damped_pendulum.cpp
    forced_pendulum.cpp
    harmonic_oscillator.cpp
    pendulum_comparison.cpp
)

add_executable(Chapter06App ${CHAPTER_SOURCES})
set_target_properties(Chapter06App PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/book_chapters)
