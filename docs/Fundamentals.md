# MML Fundamentals

**Audience:** every MML user. · **Companion docs:** [COOKBOOK.md](COOKBOOK.md) for task recipes,
subsystem docs under `docs/base|core|algorithms|systems` for depth,
[CONTRACT_POLICY.md](CONTRACT_POLICY.md) for the contributor-facing normative contract.

MML is a big library, but it is built from **five ideas that never change**. Read this page once
and every subsystem doc will feel familiar. Every code snippet below is compiled and executed by
[src/docs_demos/docs_demo_fundamentals.cpp](../src/docs_demos/docs_demo_fundamentals.cpp)
(target `MML_DocsApp`).

---

## 1. The numeric backbone: `Real`, `Complex`, and precision builds

The whole library computes in one floating-point type, chosen at compile time:

```cpp
// MMLTypeDefs.h — selected via compile definition:
//   -DMML_USE_FLOAT        → typedef float       Real;
//   -DMML_USE_DOUBLE       → typedef double      Real;   (default)
//   -DMML_USE_LONG_DOUBLE  → typedef long double Real;
typedef std::complex<Real> Complex;   // MMLBase.h — complex follows Real
```

- The single-header distribution ships all three: `MML.h` (double), `MML_float.h`,
  `MML_long_double.h`.
- Library-wide default tolerances **scale with the precision** (`MMLPrecision.h`,
  `Defaults` / `PrecisionValues<Real>`) — algorithms stay meaningful in every build.
- Write precision-portable constants as `Real(0.5)` (or the `REAL(0.5)` macro in test code),
  never bare `double` literals in generic code.

## 2. Concepts: what "a number" means (`MMLScalar` & friends)

MML uses C++20 concepts (`MMLConcepts.h`) to state *in code* which types an API accepts:

| Concept | Satisfied by | Typical use |
|---|---|---|
| `MMLReal` | `float`, `double`, `long double` | real-only algorithms |
| `MMLComplex` | `std::complex<floating-point>` | complex analysis, FFT |
| `MMLScalar` | `MMLReal` or `MMLComplex` | most numeric algorithms |
| `MMLNumeric` | any arithmetic type or complex | containers, utilities |
| `Field` | anything with `+ − * /`, `T{0}`, `T{1}` | exact/algebraic math: `Rational`, modular `ModInt`, finite fields |
| `RealFunctionCallable` | any callable `Real → Real` | passing raw lambdas to algorithms |

Because containers are templates over the element type, the *same* `Vector<T>` / `Matrix<T>`
works across numeric worlds:

```cpp
Matrix<Real>    A(3, 3);        // numerical linear algebra
Matrix<Complex> C(2, 2);        // complex eigen-analysis, FFT matrices
// and Field types (Rational, finite fields) where the math is exact

static_assert(MMLScalar<Real> && MMLScalar<Complex>);
static_assert(!MMLScalar<int>);     // int indexes things; it doesn't measure things
static_assert(Field<Complex>);
```

Concept violations fail **at compile time** with a readable requirement, not a template
back-trace.

## 3. The function-interface hierarchy

Algorithms never take raw function pointers with ad-hoc signatures — they take **interfaces**
(`mml/interfaces/IFunction.h`). This is the entry ticket to derivation, integration, root
finding, optimization, ODE solving, field operations:

| Interface | Mathematical object |
|---|---|
| `IRealFunction` | f : ℝ → ℝ |
| `IScalarFunction<N>` | f : ℝᴺ → ℝ (scalar field) |
| `IVectorFunction<N>` | f : ℝᴺ → ℝᴺ (vector field) |
| `IVectorFunctionNM<N,M>` | f : ℝᴺ → ℝᴹ |
| `IRealToVectorFunction<N>` | f : ℝ → ℝᴺ |
| `IParametricCurve<N>` | curve t ↦ ℝᴺ (extends `IRealToVectorFunction<N>`) |
| `IParametricSurface<N>` | surface (u,v) ↦ ℝᴺ |
| `I...Parametrized` variants | same, with tunable parameters (`IParametrized`) |

You supply a function in whichever form is convenient:

```cpp
// 1) Concrete adapter over a plain function (or capture-less lambda)
RealFunction f1([](Real x) { return x * x - 2; });

// 2) Adapter over std::function — lambdas WITH captures
Real shift = 2.0;
RealFunctionFromStdFunc f2(std::function<Real(Real)>(
    [shift](Real x) { return x * x - shift; }));

// 3) Or pass a raw lambda straight to an algorithm (RealFunctionCallable overloads);
//    the simple overloads return the answer directly
Real root = RootFinding::FindRootBisection([](Real x) { return x * x - 2; },
                                           0.0, 2.0, 1e-12);
```

Anything implementing the interface is a first-class citizen — an interpolated spline, a
Chebyshev approximation, or an ODE solution can be differentiated, integrated, and analyzed
exactly like a hand-written formula. **This composability is the point.**

## 4. Calling algorithms: Config in, Result out

Every iterative algorithm follows one calling pattern (reference implementation:
`mml/algorithms/RootFinding.h`):

```cpp
// Simple overload - plain answer, throws on invalid input
Real simpleRoot = RootFinding::FindRootBrent(f1, 0.0, 2.0, 1e-10);

// Full control - a Config struct
RootFinding::RootFindingConfig config;
config.tolerance      = 1e-14;
config.max_iterations = 200;
auto result = RootFinding::FindRootBrent(f1, 0.0, 2.0, config);

// Rich diagnostics — a Result struct
if (result.converged) {
    // result.root, result.function_value ≈ 0,
    // result.iterations_used, result.achieved_tolerance,
    // result.algorithm_name, result.function_evaluations
} else {
    // result.status (AlgorithmStatus enum) + result.error_message tell you WHY
}
```

The rules you can rely on:

- **Simple overloads** (trailing `xacc` tolerance) return the plain answer for throwaway use;
  the **config-based overloads** return the full Result. Both delegate to one core.
- **Config structs** (`...Config`) have sensible defaults for every field — default-construct,
  override what you need (`tolerance`, `max_iterations`, `verbose`, plus algorithm-specific
  fields).
- **Result structs** (`...Result`, extending `IterativeResultBase`) always report `converged`,
  `iterations_used`, `achieved_tolerance`, `status`, `error_message` — plus the payload
  (`root`, `eigenvalues`, `minimizer`, ...).
- **Failure to converge is data, not an exception** — check `converged`. Exceptions are
  reserved for *invalid input* (see §5).
- Some results convert implicitly to their one obvious value (`operator Real()` on
  `RootFindingResult`) — convenient, but it discards diagnostics; prefer named fields in
  anything beyond a throwaway script.

## 5. Contracts you can rely on

The user-visible guarantees, identical across the whole library (the normative, contributor
version with enforcement details is [CONTRACT_POLICY.md](CONTRACT_POLICY.md)):

| Situation | Guarantee |
|---|---|
| `a == b` | never throws; exact element-wise; different dimensions ⇒ `false` |
| tolerance comparison | named method: `IsEqualTo(b, eps)` |
| zero queries | `isZero()` is exact; `isNearZero(eps)` is tolerance-based — the name states the semantics |
| `v[i]`, `m(i,j)` | **unchecked** fast path (asserts in Debug) |
| `v.at(i)` | **checked** — throws `VectorAccessBoundsError` / `MatrixAccessBoundsError` |
| solvers/transforms | validate inputs at entry, throw MML exceptions on bad input |
| `Normalized()` | throws on a near-zero vector — never silently returns garbage; explicit `...OrZero` variants exist where a fallback is legitimate |
| any MML error | derives from `MMLException` (marker base with `message()`) *and* a matching `std::` exception — one `catch (const MMLException&)` site handles everything, and `std::exception` handlers still work |
| console output | the library **never** writes to `std::cout`; diagnostics go through streams/callbacks you provide |
| `...View` / `...Ref` names | non-owning; you guarantee the source outlives them |

```cpp
Vector<Real> v({ 1e-14, -1e-14, 0.0 });
bool exact = v.isZero();          // false — exact query
bool near  = v.isNearZero();      // true  — tolerance query (Defaults-scaled)

try {
    Real x = v.at(17);            // checked access
} catch (const MMLException& e) { // one base class catches any MML error
    // e.message() is specific and actionable
}
```

## 6. One header or many

Both styles are first-class; the namespace is always `MML`:

```cpp
#include <MML.h>                            // everything, single header (~100K lines)
// — or —
#include <mml/base/Vector/Vector.h>          // modular: pay only for what you use
#include <mml/algorithms/RootFinding.h>
```

Advanced/extension functionality (revised simplex, inferential statistics, spectral graph
analytics, heuristic optimizers, PDE solvers...) lives in the **MML-Packages** sibling library —
same namespaces, one extra include path.

---

## Where next

- [COOKBOOK.md](COOKBOOK.md) — 20 task-oriented recipes (setup → solve)
- [README_Code_examples.md](README_Code_examples.md) · [README_Usage_examples.md](README_Usage_examples.md)
- Subsystem docs: `docs/base/`, `docs/core/`, `docs/algorithms/`, `docs/systems/`
- [MIGRATION_GUIDE_1_2_to_2_0.md](MIGRATION_GUIDE_1_2_to_2_0.md) when upgrading
