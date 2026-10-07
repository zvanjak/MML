# Repository-local overlay port. For an official vcpkg registry submission,
# replace this local SOURCE_PATH with vcpkg_from_github() using a tagged release
# and the real release archive SHA512.
set(SOURCE_PATH "${CURRENT_PORT_DIR}/../..")
cmake_path(ABSOLUTE_PATH SOURCE_PATH NORMALIZE)

if(NOT EXISTS "${SOURCE_PATH}/CMakeLists.txt" OR NOT EXISTS "${SOURCE_PATH}/mml")
    message(FATAL_ERROR
        "The minimalmathlib overlay port expects ports/minimalmathlib to live inside "
        "a MinimalMathLibrary checkout. SOURCE_PATH='${SOURCE_PATH}' is not valid."
    )
endif()

vcpkg_cmake_configure(
    SOURCE_PATH "${SOURCE_PATH}"
    OPTIONS
        -DMML_BUILD_DEVELOPMENT_TARGETS=OFF
        -DMML_INSTALL=ON
)

vcpkg_cmake_install()
vcpkg_cmake_config_fixup(PACKAGE_NAME minimalmathlib CONFIG_PATH lib/cmake/minimalmathlib)

# Remove empty lib and debug directories (header-only)
file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/debug")
file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/lib")

# Install license
vcpkg_install_copyright(FILE_LIST "${SOURCE_PATH}/LICENSE.md")

# Install usage file
file(WRITE "${CURRENT_PACKAGES_DIR}/share/${PORT}/usage"
"minimalmathlib provides CMake targets:

    find_package(minimalmathlib CONFIG REQUIRED)
    target_link_libraries(your_target PRIVATE minimalmathlib::minimalmathlib)

Or simply include the single header:

    #include <MML.h>
")
