# Package Manager Support Analysis

**Issue:** MinimalMathLibrary-4e0b.22  
**Date:** 2026-01-27  
**Status:** Analysis Complete

---

## Current State

| Aspect | Status |
|--------|--------|
| **Library Type** | Header-only (single `MML.h` or modular `mml/` headers) |
| **Dependencies** | None (pure C++17) |
| **Build System** | CMake 3.20+ |
| **License** | Dual: Non-commercial free + Commercial ($50) |
| **Existing Installation** | Manual download or git clone |

---

## 1. vcpkg Support

### 1.1 What's Needed

#### A. vcpkg.json (Manifest file - in project root)

```json
{
  "name": "minimalmathlib",
  "version": "1.0.0",
  "description": "Complete C++ numerical computing toolkit - single-header library",
  "homepage": "https://github.com/zvanjak/MinimalMathLibrary",
  "license": "LicenseRef-MML-Dual",
  "supports": "!uwp"
}
```

#### B. Portfile (For vcpkg registry)

A `portfile.cmake` that:
- Downloads from GitHub release
- Copies headers to include directory
- Installs license and usage file

Location: `ports/minimalmathlib/` (in a vcpkg registry or overlay)

#### C. CMake Config File

A `minimalmathlib-config.cmake` so users can do:

```cmake
find_package(minimalmathlib CONFIG REQUIRED)
target_link_libraries(myapp PRIVATE minimalmathlib::minimalmathlib)
```

#### D. Installation Changes to CMakeLists.txt

```cmake
# Add install rules for header-only library
install(DIRECTORY mml/ DESTINATION include/mml)
install(FILES mml/single_header/MML.h DESTINATION include)
install(FILES cmake/minimalmathlib-config.cmake DESTINATION share/minimalmathlib)
```

### 1.2 vcpkg Registry Options

| Option | Effort | Reach |
|--------|--------|-------|
| **Overlay Port** | Low | Users add overlay manually |
| **Private Registry** | Medium | Can host on GitHub |
| **Official vcpkg** | High | Requires PR + review + approval |

**Recommendation:** Start with overlay port + private registry, optionally submit to official later.

---

## 2. Conan Support

### 2.1 What's Needed

#### A. conanfile.py (Recipe file)

```python
from conan import ConanFile
from conan.tools.files import copy
from conan.tools.layout import basic_layout

class MinimalMathLibraryConan(ConanFile):
    name = "minimalmathlib"
    version = "1.0.0"
    license = "LicenseRef-MML-Dual"
    url = "https://github.com/zvanjak/MinimalMathLibrary"
    description = "Complete C++ numerical computing toolkit"
    topics = ("math", "linear-algebra", "ode", "header-only")
    settings = "os", "compiler", "build_type", "arch"
    no_copy_source = True
    
    def export_sources(self):
        copy(self, "mml/*", ...)
        copy(self, "MML.h", ...)
    
    def package(self):
        copy(self, "*.h", src=..., dst="include")
    
    def package_info(self):
        self.cpp_info.bindirs = []
        self.cpp_info.libdirs = []
```

#### B. test_package/ directory

Standard Conan test package to verify installation works:

```
test_package/
├── conanfile.py
├── CMakeLists.txt
└── example.cpp
```

### 2.2 Conan Center Index Options

| Option | Effort | Reach |
|--------|--------|-------|
| **Local Recipe** | Low | Users copy conanfile.py |
| **Artifactory** | Medium | Private hosting |
| **Conan Center** | High | PR + review + CCI requirements |

**Recommendation:** Start with local recipe, optionally submit to Conan Center later.

---

## 3. Required Files Summary

| File | Purpose | Priority |
|------|---------|----------|
| `vcpkg.json` | vcpkg manifest | Required |
| `ports/minimalmathlib/portfile.cmake` | vcpkg port definition | Required |
| `ports/minimalmathlib/vcpkg.json` | vcpkg port manifest | Required |
| `cmake/minimalmathlib-config.cmake` | CMake find_package support | Required |
| `cmake/minimalmathlib-config-version.cmake` | Version checking | Optional |
| `conanfile.py` | Conan recipe | Required |
| `test_package/` | Conan test package | Required |
| `docs/INSTALLATION.md` | Installation guide | Required |
| CMakeLists.txt changes | Install rules | Required |

---

## 4. CMakeLists.txt Changes Needed

```cmake
# Add interface library target for header-only usage
add_library(minimalmathlib INTERFACE)
target_include_directories(minimalmathlib INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
)
target_compile_features(minimalmathlib INTERFACE cxx_std_17)

# Installation rules
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

install(TARGETS minimalmathlib
    EXPORT minimalmathlib-targets
)

install(DIRECTORY mml/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/mml
    FILES_MATCHING PATTERN "*.h"
)

install(FILES mml/single_header/MML.h
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(EXPORT minimalmathlib-targets
    FILE minimalmathlib-targets.cmake
    NAMESPACE minimalmathlib::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/minimalmathlib
)
```

---

## 5. License ✅

**Status:** MML now uses the **MIT License** (SPDX: `MIT`)

This is fully compatible with both vcpkg and Conan Center requirements:
- ✅ OSI-approved
- ✅ Standard SPDX identifier
- ✅ No restrictions on commercial use
- ✅ Eligible for official registry submission

---

## 6. Estimated Work Breakdown

| Task | Time |
|------|------|
| CMakeLists.txt install rules + INTERFACE target | 1h |
| cmake/minimalmathlib-config.cmake | 0.5h |
| vcpkg.json + portfile.cmake | 1.5h |
| conanfile.py + test_package | 2h |
| docs/INSTALLATION.md | 1.5h |
| Testing on Windows/Linux | 2h |
| Documentation updates (README) | 1h |
| **Total** | **~9.5h** |

---

## 7. Recommended Implementation Order

1. **CMake install rules** - Foundation for both package managers
2. **cmake/minimalmathlib-config.cmake** - Enables `find_package()`
3. **vcpkg overlay port** - Most popular on Windows
4. **conanfile.py** - Popular cross-platform
5. **docs/INSTALLATION.md** - User documentation
6. **Test on all platforms** - Verify everything works

---

## 8. Open Questions

Before implementation:

1. **Registry strategy:** Submit to official vcpkg/Conan Center eventually, or just provide overlay/local support?

2. **Version:** Use `1.0.0` or another version number?

3. **Optional dependencies:** Should the package optionally support Catch2 for building tests, or keep it strictly header-only?

4. **Scope:** Full implementation now, or start with vcpkg only?

---

## 9. File Structure After Implementation

```
MinimalMathLibrary/
├── cmake/
│   ├── minimalmathlib-config.cmake
│   └── minimalmathlib-config-version.cmake
├── ports/
│   └── minimalmathlib/
│       ├── portfile.cmake
│       └── vcpkg.json
├── test_package/
│   ├── conanfile.py
│   ├── CMakeLists.txt
│   └── example.cpp
├── conanfile.py
├── vcpkg.json
├── docs/
│   └── INSTALLATION.md
└── CMakeLists.txt (modified)
```

---

## 10. Usage Examples (Post-Implementation)

### vcpkg

```bash
# Add overlay port
vcpkg install minimalmathlib --overlay-ports=path/to/ports

# Or with manifest mode
# Add to vcpkg.json dependencies
```

```cmake
find_package(minimalmathlib CONFIG REQUIRED)
target_link_libraries(myapp PRIVATE minimalmathlib::minimalmathlib)
```

### Conan

```bash
# Create package locally
conan create . --build=missing

# Use in project
conan install . --output-folder=build
```

```cmake
find_package(minimalmathlib REQUIRED)
target_link_libraries(myapp minimalmathlib::minimalmathlib)
```

### CMake FetchContent (Alternative)

```cmake
include(FetchContent)
FetchContent_Declare(
    minimalmathlib
    GIT_REPOSITORY https://github.com/zvanjak/MinimalMathLibrary.git
    GIT_TAG v1.0.0
)
FetchContent_MakeAvailable(minimalmathlib)
target_link_libraries(myapp minimalmathlib::minimalmathlib)
```

---

## References

- [vcpkg Documentation](https://vcpkg.io/en/docs/README.html)
- [vcpkg Overlay Ports](https://vcpkg.io/en/docs/examples/overlay-triplets-linux-dynamic.html)
- [Conan 2.0 Documentation](https://docs.conan.io/2/)
- [Conan Center Index](https://github.com/conan-io/conan-center-index)
- [CMake Package Config](https://cmake.org/cmake/help/latest/manual/cmake-packages.7.html)
