# Namespaced Scalar Migration

MML's configured floating-point and complex types are now defined together in
`mml/MMLTypeDefs.h` as `MML::Real` and `MML::Complex`. The latter is exactly
`std::complex<MML::Real>`. The header includes `<complex>`; neither alias is
declared in the global namespace.

This is a clean source migration: there are no global aliases, compatibility
switches or forwarding macros. Library headers do not add global using-declarations.
An application can declare its own global `Real`, `Complex`, or `GetEnv` without
colliding with the qualified MML API.

## Consumer Code

```cpp
#include <mml/MMLTypeDefs.h>
#include <mml/base/Vector/Vector.h>

MML::Real tolerance = MML::Precision::LooseTolerance;
MML::Complex z(MML::Real(3), MML::Real(4));
MML::Vector<MML::Real> values{MML::Real(3), MML::Real(4)};
```

Inside `namespace MML` and its child namespaces, unqualified scalar names resolve
as before. Consumer code may instead import the names in its own local scope:

```cpp
{
    using MML::Real;
    using MML::Complex;
    const Real value = static_cast<Real>(1.0L / 3.0L);
    const Complex local(value, Real(0));
    std::cout << local << '\n';
}
```

Existing consumer code using `using namespace MML;` can still find the names,
but qualification or narrow local imports avoid ambiguity when the application
has identically named types.

## Literal Conversion

The global function-like `REAL(...)` macro is removed. Replace an old use with
`static_cast<MML::Real>(expression)` or `MML::Real(expression)`; within an MML/local
import scope, `Real(expression)` is sufficient. Preserve the original expression
and long-double suffix rather than rounding through a `double` literal first:

```cpp
const MML::Real precise =
    static_cast<MML::Real>(3.141592653589793238462643383279502884L);
```

The library neither defines nor consumes an application-owned `REAL` macro.
Unrelated enum members such as data-loader `ColumnType::REAL` are unchanged.

## Precision Configuration and Other Globals

`MML_USE_FLOAT`, `MML_USE_DOUBLE`, `MML_USE_LONG_DOUBLE` and CMake's `MML_REAL_TYPE`
still choose one scalar type for the whole build. Every translation unit and
linked package must use the same configuration. Namespacing is not per-translation-unit
precision selection or ABI isolation.

The former global `GetEnv` implementation is internal to
`MML::detail::visualizer`; use the public visualizer configuration APIs rather
than relying on this implementation helper. Public root/backend selection behavior
and precedence are unchanged.

Tests compile modular and standalone consumers with conflicting application-owned
global names and verify that no `REAL` macro is injected. Small syntax-only
configuration probes cover all three scalar selections without building/running
the optional full precision suites.

Runnable counterpart: [docs_demo_scalar_namespace.cpp](../src/docs_demos/base/docs_demo_scalar_namespace.cpp).
Build `MML_DocsApp` and run `--scalar-namespace`.
