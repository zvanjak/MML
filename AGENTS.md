# AI Agent Instructions for MML 2.0

This is the **public MML repository** (`zvanjak/MML`): the release-facing
distribution of Minimal Math Library.

Focus on helping users consume MML correctly, validate public examples, and keep
the public release tree coherent. Do **not** assume private development-only
folders exist here.

---

## What MML is

MML is a comprehensive C++20 numerical computing toolkit:

- single-header usage via `#include <MML.h>`
- modular-header usage via `#include <mml/...>`
- header-only core library
- no required third-party runtime dependency for the MML core
- CMake package target: `minimalmathlib::minimalmathlib`
- public package-manager path today: repository-local vcpkg overlay
- cross-platform visualizer binaries are bundled under `tools/visualizers/`

The core numerical domains include vectors, matrices, tensors, dense and sparse
linear algebra, numerical calculus, ODE/DAE solvers, root finding, optimization,
Fourier/spectral tools, geometry, differential forms, serialization, and
visualization export.

---

## Public repository boundaries

This public repo intentionally differs from the private development workspace.

Present in public MML:

- `mml/` - public headers and generated single-header variants
- `docs/` - public documentation
- `src/code_examples/` - README and documentation-backed code examples
- `src/docs_demos/` - compilable demos for documentation snippets
- `src/examples/` - larger runnable examples
- `src/visualization_examples/` - examples that use visualizers
- `tests/` - public Catch2 test suite
- `tools/visualizers/` - bundled platform visualizer binaries
- `ports/minimalmathlib/` - repository-local vcpkg overlay port
- `.github/workflows/` - public CI


---

## Include and build model

### Single-header usage

Prefer this for simple user examples:

```cpp
#include <MML.h>

int main() {
    MML::Vector<double> v{1.0, 2.0, 3.0};
    return v.size() == 3 ? 0 : 1;
}
```

Compile examples with C++20:

```bash
g++ -std=c++20 -O3 main.cpp -I/path/to/MML -o main
```

### Modular-header usage

Use modular headers when examples need targeted includes:

```cpp
#include <mml/MMLBase.h>
#include <mml/base/Vector/Vector.h>
```

Do **not** include `MML.h` and modular `mml/...` headers in the same translation
unit; the generated single header is a complete combined distribution.

### CMake usage

Installed or vcpkg/Fetched consumers should use:

```cmake
find_package(minimalmathlib CONFIG REQUIRED)
target_link_libraries(my_app PRIVATE minimalmathlib::minimalmathlib)
```

For package-only validation:

```powershell
cmake -S . -B build-package -DMML_BUILD_DEVELOPMENT_TARGETS=OFF -DMML_INSTALL=ON
cmake --build build-package --config Release --target install --parallel
```

---

## Precision model

MML's `Real` type is selected at build/configuration time:

- `MML_REAL_TYPE=float`
- `MML_REAL_TYPE=double`
- `MML_REAL_TYPE=long_double`

For ordinary local validation, use double precision:

```powershell
cmake -S . -B build -DMML_REAL_TYPE=double
```

Examples in user-facing documentation should either use explicit scalar types or
explain `Real` when precision selection matters.

---

## Build and validation commands

Routine public-release validation on Windows:

```powershell
cmake -S . -B build -DMML_REAL_TYPE=double -DMML_BUILD_BOOK=OFF -DMML_BUILD_TOOLS=OFF
cmake --build build --config Release --target MML_Tests --parallel
& .\build\tests\Release\MML_Tests.exe
```

Run the built test executable directly when validating local changes. Avoid
looping one process per test.

Optional documentation/demo validation:

```powershell
cmake --build build --config Release --target MML_DocsApp --parallel
& .\build\src\docs_demos\Release\MML_DocsApp.exe
```

Optional example aggregate:

```powershell
cmake --build build --config Release --target examples --parallel
```

---

## vcpkg overlay

The public repo currently supports a repository-local vcpkg overlay:

```powershell
git clone https://github.com/zvanjak/MML.git
vcpkg install minimalmathlib --overlay-ports=MML\ports
```

The overlay port expects `ports/minimalmathlib` to live inside a full MML
checkout. Official vcpkg registry submission is planned later and will need a
tagged release plus real source archive checksum.

---

## Visualizers

Visualizer binaries are distributed under `tools/visualizers/` for supported
platforms. Treat them as public release assets, but separate from the header-only
MML core package flow.

When changing visualizer-related docs or examples:

- keep platform distinctions clear: Windows, Linux, macOS
- do not assume a visualizer is available outside a full repository checkout
- preserve license references under `tools/visualizers/`
- validate visualization examples when the target platform and UI stack are
  available

---

## Documentation and examples

Public documentation should be runnable and consistent with the public tree.

Rules:

- Prefer snippets that use `#include <MML.h>` unless modular includes are the
  point of the example.
- Do not reference private repos, private workspace paths, or `src/book/`.
- Public README links should target `https://github.com/zvanjak/MML`.
- Substantial documentation code should have a counterpart in
  `src/docs_demos/` or `src/code_examples/`.
- If an example cannot be compiled in the public repo, fix the example or mark
  it explicitly as pseudo-code.

---

## Release hygiene

Before considering public release work done:

1. Configure and build `MML_Tests` in Release/double precision.
2. Run `build\tests\Release\MML_Tests.exe` directly.
3. Build documentation demos if documentation snippets changed.
4. Check README badges and repository links point to `zvanjak/MML`.
5. Ensure public CMake does not require `src/book/`.
6. Ensure workflows target the public branch names (`master`, `develop`).
7. Do not push to the public repository until explicitly approved by the user.

---

## Style for AI assistants

- Be precise and conservative in public-release changes.
- Preserve public API names and examples unless deliberately updating them.
- Prefer C++20 examples.
- Keep explanations user-facing; this repo is for consumers, not private
  development coordination.
- Never introduce references to private Beads IDs, private repo paths, or
  private planning documents.
