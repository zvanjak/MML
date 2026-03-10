# Installation Guide

MinimalMathLibrary (MML) is a **header-only** C++17 library. Choose your preferred installation method.

## Quick Start (Direct Include)

The simplest approach—just clone and include:

```bash
git clone https://github.com/zvonimirbandic/MinimalMathLibrary.git
```

```cpp
// Add mml/ to your include path
#include "mml/MMLBase.h"

int main() {
    MML::Vector<double, 3> v{1.0, 2.0, 3.0};
    return 0;
}
```

---

## CMake (FetchContent)

For CMake projects, use `FetchContent` to automatically download and configure:

```cmake
cmake_minimum_required(VERSION 3.20)
project(MyProject)

include(FetchContent)

FetchContent_Declare(
    minimalmathlib
    GIT_REPOSITORY https://github.com/zvonimirbandic/MinimalMathLibrary.git
    GIT_TAG        main  # or specific version tag
)

FetchContent_MakeAvailable(minimalmathlib)

add_executable(my_app main.cpp)
target_link_libraries(my_app PRIVATE MML::minimalmathlib)
```

### CMake Install (System-Wide)

Build and install MML system-wide:

```bash
git clone https://github.com/zvonimirbandic/MinimalMathLibrary.git
cd MinimalMathLibrary
cmake -B build -DCMAKE_INSTALL_PREFIX=/usr/local
cmake --install build
```

Then use `find_package`:

```cmake
find_package(minimalmathlib REQUIRED)
target_link_libraries(my_app PRIVATE MML::minimalmathlib)
```

---

## vcpkg

### Manifest Mode (Recommended)

Add to your project's `vcpkg.json`:

```json
{
    "dependencies": ["minimalmathlib"]
}
```

If MML isn't in the vcpkg registry yet, use the overlay port from this repository:

```bash
# Clone MML for its overlay port
git clone https://github.com/zvonimirbandic/MinimalMathLibrary.git

# Configure with overlay
cmake -B build -S . \
    -DCMAKE_TOOLCHAIN_FILE=[vcpkg-root]/scripts/buildsystems/vcpkg.cmake \
    -DVCPKG_OVERLAY_PORTS=path/to/MinimalMathLibrary/ports
```

### Classic Mode

```bash
# Using overlay port
vcpkg install minimalmathlib --overlay-ports=path/to/MinimalMathLibrary/ports

# Once in registry
vcpkg install minimalmathlib
```

### CMake Integration

```cmake
find_package(minimalmathlib CONFIG REQUIRED)
target_link_libraries(my_app PRIVATE MML::minimalmathlib)
```

---

## Conan

### Add Dependency

**conanfile.txt:**
```ini
[requires]
minimalmathlib/1.0.0

[generators]
CMakeDeps
CMakeToolchain
```

**conanfile.py:**
```python
from conan import ConanFile
from conan.tools.cmake import cmake_layout

class MyProjectConan(ConanFile):
    settings = "os", "compiler", "build_type", "arch"
    generators = "CMakeDeps", "CMakeToolchain"
    
    def requirements(self):
        self.requires("minimalmathlib/1.0.0")
    
    def layout(self):
        cmake_layout(self)
```

### Install Dependencies

```bash
conan install . --output-folder=build --build=missing
```

### CMake Integration

```cmake
cmake_minimum_required(VERSION 3.20)
project(MyProject)

find_package(minimalmathlib REQUIRED)

add_executable(my_app main.cpp)
target_link_libraries(my_app PRIVATE MML::minimalmathlib)
```

### Build

```bash
cd build
cmake .. -DCMAKE_TOOLCHAIN_FILE=conan_toolchain.cmake
cmake --build .
```

### Local Export (Development)

To use MML from local source before it's in Conan Center:

```bash
cd MinimalMathLibrary
conan export .
conan install --requires=minimalmathlib/1.0.0 --build=missing
```

---

## Verification

Test your installation with this example:

```cpp
#include <iostream>
#include "mml/MMLBase.h"

int main() {
    // Vector operations
    MML::Vector<double, 3> v1{1.0, 2.0, 3.0};
    MML::Vector<double, 3> v2{4.0, 5.0, 6.0};
    double dot = v1.dot(v2);
    std::cout << "Dot product: " << dot << std::endl;
    
    // Matrix operations
    MML::Matrix<double, 2, 2> A{{1, 2}, {3, 4}};
    MML::Vector<double, 2> x{1.0, 2.0};
    auto b = A * x;
    std::cout << "Matrix-vector product: [" << b[0] << ", " << b[1] << "]" << std::endl;
    
    // Linear solve
    MML::Matrix<double, 3, 3> M{{4, 12, -16}, {12, 37, -43}, {-16, -43, 98}};
    MML::Vector<double, 3> rhs{1.0, 2.0, 3.0};
    
    MML::LUSolver<double, 3> solver;
    auto solution = solver.solve(M, rhs);
    std::cout << "Solution: [" << solution[0] << ", " << solution[1] << ", " << solution[2] << "]" << std::endl;
    
    return 0;
}
```

---

## Requirements

- **C++17** or later
- **CMake 3.20+** (for CMake-based installation)
- **vcpkg** (for vcpkg installation)
- **Conan 2.0+** (for Conan installation)

---

## Troubleshooting

### "Cannot find minimalmathlib"

- **vcpkg**: Ensure overlay port path is correct and `CMAKE_TOOLCHAIN_FILE` is set
- **Conan**: Run `conan export .` in MML directory first
- **CMake**: Verify `CMAKE_PREFIX_PATH` includes MML install location

### Compiler Errors

- Ensure C++17 is enabled: `-std=c++17` or `/std:c++17`
- Check that all include paths are correctly set

### Windows-Specific

- Use forward slashes in paths or escape backslashes
- Consider using CMake presets for consistent configuration

---

## Package Manager Status

| Method | Status | Notes |
|--------|--------|-------|
| Direct Include | ✅ Ready | Just clone and use |
| CMake FetchContent | ✅ Ready | Recommended for CMake projects |
| CMake Install | ✅ Ready | System-wide installation |
| vcpkg (overlay) | ✅ Ready | Use overlay port from this repo |
| vcpkg (registry) | 🔄 Pending | Submission to vcpkg planned |
| Conan (local) | ✅ Ready | Use `conan export .` |
| Conan Center | 🔄 Pending | Submission planned |

---

## License

MinimalMathLibrary is released under the [MIT License](../LICENSE.md).
