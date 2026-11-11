# Standard Functions - Mathematical Function Wrappers

**File**: `mml/base/StandardFunctions.h`

Unified namespace for standard mathematical functions.

MML's C++17-and-later special-function wrappers do not require
`__STDCPP_WANT_MATH_SPEC_FUNCS__`. Modular and generated headers do not define
that legacy opt-in macro, and including `<cmath>` before MML is supported.
Automatic selection reads `__cpp_lib_math_special_functions` after including
`<cmath>`: a value of at least `201603L` selects the native wrappers; an absent or
older feature value selects the existing portable fallback subset. Compiler
version and operating system do not imply availability. In particular, Clang
with libstdc++ can select native functions while the same compiler with libc++
can select the fallback subset, including on Linux.

For explicit portability testing, define `MML_HAS_STD_SPECIAL_FUNCTIONS=0` before
including MML to force its existing fallback subset. Keep this selection
consistent within a program; it does not create substitutes for every native-only
wrapper listed on this page. An explicit `MML_HAS_STD_SPECIAL_FUNCTIONS=1` also
takes precedence, but requires a library that actually supplies the native
functions; it cannot add missing library facilities. The separate `MML_StandardFallback_Tests` and
`MML_SingleHeader_StandardFallback_Tests` targets force this path for modular and
generated double headers, independently of the normal test executables.
Complete elliptic E (`Comp_ellint_2`) always uses MML's shared AGM implementation,
regardless of this selection, to avoid native libstdc++ accuracy loss near unit modulus.

## Overview

`MML::Functions` namespace provides consistent access to common scalar mathematical functions, wrapping standard library implementations.

**Key Features:**
- ✅ **Consistent naming** - Capitalized function names (Sin, Cos, Exp, etc.)
- ✅ **Type-safe wrappers** - All functions operate on `Real` type
- ✅ **Special functions** - Gamma, error, Bessel, Legendre, etc.
- ✅ **Inline implementations** - Zero overhead

---

## Trigonometric Functions

```cpp
using namespace MML::Functions;

// Basic trig
Real s = Sin(Constants::PI / 6);    // 0.5
Real c = Cos(Constants::PI / 3);    // 0.5
Real t = Tan(Constants::PI / 4);    // 1.0

// Inverse trig
Real angle = Asin(0.5);             // π/6
Real angle2 = Acos(0.5);            // π/3
Real angle3 = Atan(1.0);            // π/4
```

---

## Hyperbolic Functions

```cpp
// Hyperbolic
Real sh = Sinh(1.0);    // (e - 1/e) / 2
Real ch = Cosh(1.0);    // (e + 1/e) / 2
Real th = Tanh(1.0);    // sh / ch

// Inverse hyperbolic
Real x1 = Asinh(1.0);
Real x2 = Acosh(2.0);
Real x3 = Atanh(0.5);
```

---

## Exponential and Logarithmic

```cpp
// Exponential
Real e_x = Exp(1.0);      // e ≈ 2.71828
Real pow_val = Pow(2.0, 10.0);  // 1024

// Logarithmic
Real ln = Log(Constants::E);    // 1.0 (natural log)
Real log10_val = Log10(100.0);  // 2.0
```

---

## Roots and Powers

```cpp
// Square root
Real sqrt_val = Sqrt(16.0);      // 4.0

// Power (use for nth root)
Real pow_val = Pow(2.0, 0.5);    // sqrt(2)
Real cbrt_val = Pow(27.0, 1.0/3.0);  // 3.0 (cube root)
```

---

## Special Functions

### Error Functions
```cpp
// Error function: erf(x) = (2/√π) ∫₀ˣ e^(-t²) dt
Real erf_val = Erf(1.0);         // ≈ 0.8427

// Complementary error function: erfc(x) = 1 - erf(x)
Real erfc_val = Erfc(1.0);       // ≈ 0.1573
```

### Gamma Functions
```cpp
// Gamma function: Γ(x) = ∫₀^∞ t^(x-1) e^(-t) dt
Real gamma = TGamma(5.0);        // 4! = 24

// Log-gamma: ln(Γ(x))
Real lgamma = LGamma(5.0);       // ln(24)
```

#### Portable incomplete gamma functions

The separate `<mml/base/SpecialFunctions.h>` header provides
`MML::SpecialFunctions::IncompleteGammaLower(a, x)` and `IncompleteGammaUpper(a, x)`.
They evaluate the defining lower and upper integrals, not probabilities.
Log-domain scaling avoids multiplying an underflowed regularized probability by
an overflowing complete gamma. For example, both results below are finite even
though `Gamma(200)` overflows double precision:

```cpp
#include <mml/base/SpecialFunctions.h>
#include <iostream>

void Docs_Demo_IncompleteGamma()
{
    using namespace MML;
    const Real lower = SpecialFunctions::IncompleteGammaLower(200, 1);
    const Real upper = SpecialFunctions::IncompleteGammaUpper(200, 2000);
    std::cout << lower << '\n';  // approximately 0.001848593963122709
    std::cout << upper << '\n';  // approximately 2.298767256466363e-212
}
```

Both functions require finite `a > 0` and `x >= 0` (including positive infinity).
At `x = 0` they return `0` and `Gamma(a)` respectively; at positive infinity
they return `Gamma(a)` and `0`. Genuine final overflow returns positive infinity;
underflow rounds to a subnormal value or zero, following the active floating-point
environment. Neither endpoint forms `0 * infinity`.

The implementation shares the regularized gamma family's lower series
(`x < a + 1`) and upper continued fraction. A complementary result uses `expm1`
in log space. Invalid domains throw `DomainError`; exhausted
`Precision::IncompleteGammaMaxIterations` throws `ConvergenceError`;
arithmetic failures throw `NumericalMethodError`. Large shapes use the bounded
density integration described below instead of exhausting the lower series.
For tiny shapes and `x <= 1`, an alternating lower-gamma series factors out
the shape before taking the complement, retaining finite upper incomplete
gamma values even when the regularized upper tail rounds to zero.

Runnable counterpart:
[`Docs_Demo_IncompleteGamma`](../../src/docs_demos/base/docs_demo_standard_functions.cpp),
available through `MML_DocsApp --gamma-scaling`.

#### Stable forward gamma/beta normalization and large shapes

The normalization and density-integration paths use at least double precision internally.
The gamma series and continued fraction retain that working precision, including
their convergence comparisons; conversion to `Real` occurs at result and diagnostic boundaries.
Stirling corrections and a cancellation-free `log(1+r)-r` series evaluate
gamma deviance near the mode. Positive-shape `Beta`/`LnBeta` use stable
gamma ratios or Stirling normalization; small-shape normalization uses
`log(Gamma(1+a))`, without first rounding `1+a` to one. The incomplete-beta
fraction forms ratios before products to avoid artificial overflow.

At `a >= Precision::ForwardShapeQuadratureSwitch` (1,000), gamma tails
integrate the exact density relative to its endpoint value. Beta uses the
same approach when both shapes are at least
`Precision::ForwardShapeStirlingMin` (16) and their sum reaches the switch.
The integration coordinate is scaled by the smaller of the local density
width and inverse log-density slope, so the interval does not grow with
shape. An adaptive eight-point Gauss rule compares each interval with its
two halves. The logarithmic tail window is
`-Precision::ForwardShapeQuadratureTailLogs * log(epsilon)` in this coordinate.
This is density integration, not a normal approximation substituted for the CDF.

`Precision::ForwardShapeQuadratureMaxEvaluations` (50,000) and
`Precision::ForwardShapeQuadratureMaxDepth` (20) bound the work.
Exhaustion throws `ConvergenceError`; invalid arithmetic throws
`NumericalMethodError`. The legacy series/fraction budgets remain enforced
when those methods are selected. Incomplete-beta interiors still reject an
unrepresentable shape sum explicitly.

`IncompleteBeta` combines the log probability and log normalization before
exponentiating, avoiding an underflowed probability multiplied by an
overflowing complete beta. Its endpoints are zero and `Beta(a,b)`.
Genuine final overflow/underflow follows the floating-point range.

```cpp
const Real gammaP = SpecialFunctions::RegularizedGammaP(Real(1e6), Real(1e6));
const Real betaP = SpecialFunctions::RegularizedBetaI(Real(0.5), Real(1e12), Real(1e12));
const Real gammaMedian = SpecialFunctions::RegularizedGammaPInv(Real(1e6), Real(0.5));
const Real finiteBeta = SpecialFunctions::IncompleteBeta(
    Real(0.5), 2, std::numeric_limits<Real>::denorm_min());
std::cout << gammaP << " " << betaP << " " << gammaMedian << " " << finiteBeta << '\n';
```

In double precision these are approximately `0.500132980760873`, `0.5`,
`999999.666666686`, and `0.193147180559945`. Independent high-precision
density-integral regressions cover central shapes through `1e12`, disparate
and tiny shapes. This is measured coverage, not a uniform accuracy guarantee
over every finite parameter. Quantile ULP spacing can dominate probability
roundtrips at very large shapes; unresolved numerical cases still fail explicitly.

Runnable counterpart:
[`Docs_Demo_ForwardShapes`](../../src/docs_demos/base/docs_demo_standard_functions.cpp),
available through `MML_DocsApp --forward-shapes`.

#### Portable inverse lower gamma

`SpecialFunctions::RegularizedGammaPInv(a, p)` solves `P(a, x) = p`.
Its domain is **finite `a > 0` and finite `p` in `[0,1]`**. Invalid inputs,
including NaN, infinite shapes/probabilities and out-of-range probabilities,
throw `DomainError`; probabilities are no longer silently clamped.
The shape is validated even at the exact probability endpoints.

For valid shapes, `p = 0` returns positive zero and `p = 1` returns positive
infinity. Positive quantiles below the representable range round to a subnormal
or zero, rather than an artificial minimum iterate.

```cpp
const Real lower = SpecialFunctions::RegularizedGammaPInv(2, Real(1e-100));
const Real tiny = SpecialFunctions::RegularizedGammaPInv(Real(0.5), Real(1e-200));
const Real p = std::nextafter(Real(1), Real(0));
const Real upper = SpecialFunctions::RegularizedGammaPInv(1, p);
std::cout << std::setprecision(std::numeric_limits<Real>::max_digits10);
std::cout << lower << '\n' << tiny << '\n' << upper << '\n';
```

In double precision the results are approximately `1.414213562373e-50`,
`0`, and `36.7368005696771`. Shapes `1` and `1/2` use the exact exponential
and error-function identities. Other shapes use a positive bracket and
safeguarded Newton steps against `log(P)` or `log(Q)`, selecting the upper
tail for `p > 1/2`. Near zero, the power-law limit is evaluated in log space;
small-shape `log(Gamma(1+a))` uses its zeta-coefficient Taylor expansion to
avoid losing `a` in `1+a`.

The leading lower-tail quantile is returned directly when its first relative
correction, `x/(a+1)`, is at most `Precision::Eps`. This avoids making an already
accurate tiny quantile worse through coarsely spaced log-tail iterates, particularly
with `Real=float`.

`Precision::IncompleteGammaInverseMaxIterations` bounds log-tail evaluations,
including bracket expansion. Outside that asymptotic region, both a relative quantile bracket and a log-tail
residual check are required. The latter includes kernel rounding and the
probability change across a quantile ULP; it is not a uniform machine-relative
probability guarantee. Exhaustion or unresolved stagnation throws
`ConvergenceError`, never an unchecked final iterate. Forward gamma
convergence/arithmetic failures propagate. The inverse shares the stable
forward normalization and large-shape integration path above, including the
previously failing `a = 1e6, p = 1/2` case.

Runnable counterpart:
[`Docs_Demo_GammaInverseTails`](../../src/docs_demos/base/docs_demo_standard_functions.cpp),
available through `MML_DocsApp --gamma-inverse-tails`. The demo also verifies
strict rejection of an out-of-range probability.

#### Portable inverse incomplete beta

`SpecialFunctions::RegularizedBetaIInv(p, a, b)` solves `I_x(a,b) = p`.
The argument order is **probability first**. Both shapes must be finite and
strictly positive, and `p` must be finite and in `[0,1]`. Invalid inputs throw
`DomainError`, including invalid shapes at the exact probability endpoints.
Out-of-range probabilities are no longer silently clamped.

For valid shapes, `p = 0` returns positive zero and `p = 1` returns one.
Mathematically interior quantiles can round to subnormals, zero, or one;
there is no artificial interior floor.

```cpp
const Real lower = SpecialFunctions::RegularizedBetaIInv(Real(1e-100), 2, 1);
const Real uniform = SpecialFunctions::RegularizedBetaIInv(Real(1e-300), 1, 1);
const Real tiny = SpecialFunctions::RegularizedBetaIInv(Real(1e-200), Real(0.5), Real(0.5));
const Real rounded = SpecialFunctions::RegularizedBetaIInv(Real(0.99), 2, Real(0.1));
std::cout << std::setprecision(std::numeric_limits<Real>::max_digits10);
std::cout << lower << '\n' << uniform << '\n' << tiny << '\n' << rounded << '\n';
```

In double precision these results are approximately `1e-50`, `1e-300`, `0`,
and `1`. Uniform and unit-shape distributions use their closed-form identities,
with endpoint tie corrections to avoid premature zero/one from intermediate rounding.
General shapes use safeguarded Newton steps against the smaller probability
tail, in the log distance from the nearer support boundary. Endpoint rounding
is checked at half an endpoint ULP rather than inferred from an initial guess.

For a normal, non-reflected lower-tail quantile, the leading power-law result
is accepted when its first relative correction, `abs(b-1)*x/(a+1)`, is at most
`Precision::Eps`. Subnormal quantiles and reflected endpoint rounding retain
the checked bracket path.

`Precision::IncompleteBetaInverseMaxIterations` bounds all log-tail evaluations.
Outside the asymptotic shortcut, both a relative endpoint-distance bracket and a scale-aware log-tail residual
are checked. Exhaustion or unresolved stagnation throws `ConvergenceError`,
never an unchecked final iterate. The inverse shares the stable normalization,
ratio-ordered continued fraction, and large-shape density integration above.
Their convergence and arithmetic failures propagate. Shapes `a = b = 1e12`
now resolve the median instead of exhausting the legacy fraction.

Accuracy and representability follow the configured `Real` type. The tiny
decimal inputs in the examples above are double-precision examples; an input
that has already rounded to zero in float is the exact probability endpoint,
not a representable tail. Wider internal arithmetic cannot recover precision
already lost in an input parameter or in the returned quantile.
Inverse roundtrips must account for the documented log-tail residual and the
probability change across a quantile ULP, rather than assuming a uniform relative
probability tolerance. Similarly, rounding `erfc(x)` near one can dominate a
small `ErfcInv` roundtrip through its inverse derivative.
Inverse-beta symmetry comparisons must propagate the log-tail error through the inverse
density for both solves in every precision, including `double` and `long double`.

This portable inverse's endpoint-rounding contract is distinct from
`Statistics::BetaDistribution::inverseCdf`, which retains its checked,
resolvable-interior quantile contract and existing probability clamping.

Runnable counterpart:
[`Docs_Demo_BetaInverseTails`](../../src/docs_demos/base/docs_demo_standard_functions.cpp),
available through `MML_DocsApp --beta-inverse-tails`. The demo also checks
strict rejection of an out-of-range probability.

#### Portable real Lambert W branches

The same portable header provides `SpecialFunctions::LambertW0(x)` and
`SpecialFunctions::LambertWm1(x)`, solving `w * exp(w) = x`.
The principal branch is increasing on `[-1/e, +infinity)` with `w >= -1`;
the secondary branch is decreasing on `[-1/e, 0)` with `w <= -1`.

Only the rounded `Real` representation of `-1/e` is treated as the exact
endpoint and returns `-1` on both branches. Larger valid inputs are not
flattened by a tolerance. In particular, the first representable neighbor
above that endpoint already has two distinct real results:

```cpp
#include <mml/base/SpecialFunctions.h>
#include <iostream>
#include <iomanip>

void Docs_Demo_LambertW()
{
    using namespace MML;
    const Real endpoint = Real(-0.367879441171442321595523770161460867L);
    const Real x = std::nextafter(endpoint, Real(0));
    const Real principal = SpecialFunctions::LambertW0(x);
    const Real secondary = SpecialFunctions::LambertWm1(x);
    std::cout << std::setprecision(std::numeric_limits<Real>::max_digits10);
    std::cout << principal << '\n' << secondary << '\n';
}
```

In double precision these results are approximately `-0.9999999846957459`
and `-1.0000000153042543`. At the original regression input
`endpoint + 5e-11`, they are approximately `-0.9999835128792688` and
`-1.0000164873019500`. The exact binary input matters: this inverse is
ill-conditioned near the branch point, and a small defining-equation residual
alone cannot establish accurate branch values.

The shared implementation uses a split `-1/e` constant and a ninth-order
Puiseux expansion near the endpoint, following the branch expansion in
Corless et al., *On the Lambert W Function* (1996).
Outside that small neighborhood it refines the logarithmic equation
`w + log(abs(w)) = log(abs(x))`, avoiding exponential overflow and underflow.
The expansion's `Precision::LambertWSeriesLimit` is a method-selection limit,
not a tolerance that changes the domain.

Both functions reject NaN and inputs below the rounded endpoint with
`DomainError`. The principal branch returns zero at either signed zero and
positive infinity at positive infinity. The secondary branch rejects either
signed zero, positive inputs and both infinities. Arithmetic or branch-loss
failures throw `NumericalMethodError`; exhausted
`Precision::LambertWMaxIterations` throws `ConvergenceError` instead of returning
an unchecked iterate.

Runnable counterpart:
[`Docs_Demo_LambertW`](../../src/docs_demos/base/docs_demo_standard_functions.cpp),
available through `MML_DocsApp --lambertw`.

#### Dawson and Ei at large arguments

`SpecialFunctions::Dawson(x)` is defined for every finite real argument and is
odd. Near zero it retains the Maclaurin approximation; at ordinary magnitudes
it retains the six-term Numerical Recipes sampling approximation (approximately
`Precision::LooseTolerance` relative accuracy). For
`abs(x) >= Precision::DawsonAsymptoticLimit`, a reciprocal asymptotic expansion
is evaluated **before** any integer grid index is formed. This avoids
out-of-range integer conversions, overflowing `x*x` or `2*x`, and preserves
representable subnormal tails.

`SpecialFunctions::ExpIntegralEi(x)` retains its nonzero real domain and
`Ei(x) = -E1(-x)` for negative inputs. Its positive large-argument expansion
uses split exponential scaling rather than forming `exp(x)` first.
For the default double precision, `Ei(710)` is finite even though `exp(710)`
overflows:

```cpp
#include <mml/base/SpecialFunctions.h>
#include <iostream>
#include <iomanip>

void Docs_Demo_SpecialFunctionRanges()
{
    using namespace MML;
    const Real tail = SpecialFunctions::Dawson(std::numeric_limits<Real>::max());
    const Real ei = SpecialFunctions::ExpIntegralEi(710);
    std::cout << std::setprecision(std::numeric_limits<Real>::max_digits10);
    std::cout << tail << '\n' << ei << '\n';
}
```

In double precision the results are approximately `2.781342323134e-309`
and `3.150915688206201e305`. The maximum finite positive Ei input is about
`716.3554905424517`; the next representable input already has an overflowing
result. Genuine final overflow returns positive infinity. Underflow follows
the active floating-point environment, rounding to subnormals or signed zero.
The old sampling approximation is not made uniformly machine-accurate:
the method transitions retain its ordinary approximation accuracy.

Both functions reject NaN with `DomainError`. Dawson preserves either signed
zero and returns signed zero at the corresponding positive/negative infinity.
Ei still rejects either signed zero; its positive-infinity limit is positive
infinity and its negative-infinity limit is negative zero. The positive Ei
series and asymptotic expansion have the named
`Precision::ExpIntegralEiMaxIterations` budget; the latter uses least-term
truncation when its terms begin growing. Dawson's tail uses
`Precision::DawsonMaxTerms` and reports `ConvergenceError` if it cannot meet
its relative stopping criterion. Exhausted positive Ei budgets also throw
`ConvergenceError`, rather than silently returning an unchecked partial sum.
The separate negative-input `E1`/`En` kernels are unchanged by this repair.

Runnable counterpart:
[`Docs_Demo_SpecialFunctionRanges`](../../src/docs_demos/base/docs_demo_standard_functions.cpp),
available through `MML_DocsApp --special-function-ranges`.

### Sine and cosine integrals

`SpecialFunctions::SinCosIntegral(x, si, ci)` computes both integrals;
`SineIntegral(x)` and `CosineIntegral(x)` use the same kernel. For finite
negative inputs MML keeps a **real-valued** extension: Si is odd, while
Ci is even and uses `log(abs(x))`. This is not the complex principal branch
of Ci on the negative real axis.

| Input | Si | Ci |
|-------|----|----|
| `+0` / `-0` | Same signed zero | Negative infinity |
| Positive / negative infinity | Positive / negative `pi/2` | Positive zero |
| NaN | `DomainError` | `DomainError` |

The paired call rejects NaN before modifying either output. Only exact
zero takes the singular endpoint branch: nonzero finite values, including
subnormals, retain their finite logarithmic Ci value.

```cpp
const Real infinity = std::numeric_limits<Real>::infinity();
Real si, ci;
SpecialFunctions::SinCosIntegral(-Real(0), si, ci);
std::cout << std::signbit(si) << ' ' << ci << '\n';  // 1 -inf
const Real sineLimit = SpecialFunctions::SineIntegral(-infinity);
const Real cosineLimit = SpecialFunctions::CosineIntegral(infinity);
std::cout << sineLimit << ' ' << cosineLimit << '\n';  // -pi/2 0
```

The finite-input power series (`abs(x) <= 2`), continued fraction, method
cutoffs and iteration limits are unchanged by the endpoint repair; this
does not establish a new finite-range accuracy or convergence guarantee.

Runnable counterpart:
[`Docs_Demo_SineCosineEndpoints`](../../src/docs_demos/base/docs_demo_standard_functions.cpp),
available through `MML_DocsApp --sici-endpoints`. The demo also checks
NaN rejection and preservation of the paired outputs.

### Bessel Functions
```cpp
// Spherical Bessel function of first kind
Real j_0 = SphBessel(0, 1.0);
Real j_1 = SphBessel(1, 1.0);
```

### Legendre Polynomials
```cpp
// Legendre polynomial P_n(x)
Real P_0 = Legendre(0, 0.5);     // 1
Real P_1 = Legendre(1, 0.5);     // 0.5
Real P_2 = Legendre(2, 0.5);     // -0.125

// Spherical harmonics (associated Legendre)
Real Y = SphLegendre(2, 1, 0.5);
```

### Hermite and Laguerre Polynomials
```cpp
// Hermite polynomial H_n(x)
Real H_0 = Hermite(0, 1.0);      // 1
Real H_1 = Hermite(1, 1.0);      // 2
Real H_2 = Hermite(2, 1.0);      // 2

// Laguerre polynomial L_n(x)
Real L_0 = Laguerre(0, 1.0);     // 1
Real L_1 = Laguerre(1, 1.0);     // 0
Real L_2 = Laguerre(2, 1.0);     // -0.5
```

### Elliptic Integrals
```cpp
// Complete elliptic integral of first kind
Real K = Comp_ellint_1(0.5);

// Complete elliptic integral of second kind
Real E = Comp_ellint_2(0.5);
```

### Portable Legendre and elliptic-E conventions

`AssocLegendre(l,m,x)` uses the standard-library convention **without** the
Condon-Shortley factor `(-1)^m`. `SphLegendre(l,m,theta)` retains that factor
and the spherical-harmonic normalization. The associated polynomial accepts
`x` in `[-1,1]`; for a valid finite argument, `m > l` returns zero. Its NaN
argument propagates, and `|x| > 1` (including infinities) throws MML `DomainError`
in both native and fallback selection.

`Comp_ellint_2(k)` takes the **modulus** `k`, not the parameter `k*k`. It is even,
has `E(0)=pi/2` and exact endpoint values `E(1)=E(-1)=1`. NaN propagates and
`|k| > 1` (including infinities) throws MML `DomainError` in either selection.
Both selections use the same AGM correction weights `k*k/2, c1*c1, 2*c2*c2, ...`,
a relative AGM gap criterion, and the type-scaled
`Precision::EllipticAGMMaxIterations` budget. Exhaustion throws `ConvergenceError`.
The first-kind elliptic integral and other fallback algorithms are not rewritten.

The following example is checked by `Docs_Demo_StandardFallbackParity()` in
`src/docs_demos/base/docs_demo_standard_fallback.cpp`, runnable through
`MML_DocsApp --standard-fallback-parity`. Its same source is also compiled and
executed under forced fallback during validation:

```cpp
const Real associated = Functions::AssocLegendre(1, 1, 0);
const Real spherical = Functions::SphLegendre(1, 1, Constants::PI / 2);
const Real elliptic = Functions::Comp_ellint_2(Real(0.5));
std::cout << associated << '\n';  // 1: no Condon-Shortley phase
std::cout << spherical << '\n';   // approximately -0.345494
std::cout << elliptic << '\n';    // approximately 1.467462209339427
```

The parity regressions cover analytic low-degree values, native-library
comparisons through degree 12, independent complete-E values and endpoint neighbors
under both selections, and the explicit
domain policy. They do not establish uniform accuracy at arbitrary degree or
for every tiny spherical-harmonic value near an angular pole.

### Riemann Zeta Function
```cpp
// ζ(s) = Σ(n=1 to ∞) 1/n^s
Real zeta_2 = RiemannZeta(2.0);  // π²/6 ≈ 1.6449
```

---

## Complete Function List

### Trigonometric
- `Sin(x)`, `Cos(x)`, `Tan(x)`
- `Sec(x)`, `Csc(x)`, `Ctg(x)` (reciprocals)
- `Asin(x)`, `Acos(x)`, `Atan(x)`

### Hyperbolic
- `Sinh(x)`, `Cosh(x)`, `Tanh(x)`
- `Sech(x)`, `Csch(x)`, `Ctgh(x)` (reciprocals)
- `Asinh(x)`, `Acosh(x)`, `Atanh(x)`

### Exponential/Logarithmic
- `Exp(x)`, `Pow(x, y)`
- `Log(x)`, `Log10(x)`

### Roots
- `Sqrt(x)`
- Use `Pow(x, 1.0/n)` for nth root

### Special Functions
- `Erf(x)`, `Erfc(x)` - Error functions
- `TGamma(x)`, `LGamma(x)` - Gamma functions
- `SphBessel(n, x)` - Spherical Bessel
- `CylBesselJ(n, x)`, `CylBesselI(n, x)`, `CylBesselK(n, x)` - Cylindrical Bessel
- `CylBesselY(n, x)` / `CylNeumann(n, x)` - Neumann functions
- `Legendre(n, x)` - Legendre polynomials
- `AssocLegendre(l, m, x)` - Associated Legendre
- `SphLegendre(l, m, x)` - Spherical harmonics
- `Hermite(n, x)` - Hermite polynomials
- `Laguerre(n, x)`, `AssocLaguerre(n, m, x)` - Laguerre polynomials
- `Comp_ellint_1(k)`, `Comp_ellint_2(k)`, `Comp_ellint_3(k, n)` - Complete elliptic integrals
- `Ellint_1(k, phi)`, `Ellint_2(k, phi)`, `Ellint_3(k, n, phi)` - Incomplete elliptic integrals
- `RiemannZeta(s)` - Riemann zeta function
- `Beta(x, y)` - Beta function
- `Expint(x)` - Exponential integral

### Complex reciprocal wrappers

All six reciprocals (`Sec`, `Csc`, `Ctg`, `Sech`, `Csch`, `Ctgh`) accept
`Complex` and `std::complex<float/double/long double>`, as well as the corresponding
real scalar types. Their return type remains the input type.

`Sec`, `Csc`, `Ctg`, `Csch`, and `Ctgh` throw `DomainError` when the magnitude of
their computed denominator is **strictly less than** machine epsilon of its real
magnitude type. For a complex input this is one complex magnitude, not a separate
test of each component; equality with epsilon is not rejected. This retains the
existing real-input guard rather than introducing a different singularity cutoff.
`Sech` retains its unguarded `1/cosh(z)` behavior, including near complex poles.
Non-finite inputs and results continue to follow the underlying standard-library
operations; this repair adds no separate NaN or infinity policy.

This example is compiled and checked by `Docs_Demo_ComplexReciprocals()` in
`src/docs_demos/base/docs_demo_standard_functions.cpp`, runnable through
`MML_DocsApp --complex-reciprocals`:

```cpp
const Complex z(1, 1);
const Complex values[] = {
    Functions::Sec(z), Functions::Csc(z), Functions::Ctg(z),
    Functions::Sech(z), Functions::Csch(z), Functions::Ctgh(z)
};
for (const Complex value : values)
    std::cout << value << '\n';

bool rejected = false;
try {
    Functions::Csc(Complex(0));
}
catch (const DomainError&) {
    rejected = true;
    std::cout << "Complex cosecant rejected its zero denominator\n";
}
```

### Utility
- `Factorial(n)` - n! as Real; `DomainError` for negative input, positive infinity on overflow
- `FactorialInt(n)` - exact n! as long long; `DomainError` for negative input or an unrepresentable result
- `FactorialStirling(n)` - Stirling's approximation

### Floating factorial range

`Factorial` retains the ascending multiplication sequence for finite results.
It checks the active `Real` range before each multiplication and returns positive
infinity as soon as the result would overflow, without evaluating the remaining
factors. This also holds under directed rounding; overflow is not returned as a
saturated finite value. Even `INT_MAX` terminates at the range boundary, rather
than running a loop whose signed counter could overflow. With double precision,
170! is finite and 171! returns positive infinity; the cutoff is derived from
the actual `Real` range, not a hard-coded double-precision limit.

This example is covered by `Docs_Demo_IntegerBoundaries()` and is runnable through
`MML_DocsApp --integer-boundaries`:

```cpp
const Real ordinary = Functions::Factorial(20);
const Real huge = Functions::Factorial(std::numeric_limits<int>::max());
std::cout << ordinary << '\n';
std::cout << huge << '\n';  // positive infinity
```

### Exact integer boundaries

`FactorialInt` uses the shared prechecked exact-integer multiplication kernel.
On the supported signed 64-bit representation, 0! through 20! are exact; every
`n >= 21`, including `INT_MAX`, throws `DomainError` before a multiplication can
overflow. No wrapped, saturated, or floating-point substitute is returned.
The separate floating-point `Factorial` returns positive infinity on overflow,
not an exact-integer substitute. `FactorialStirling` remains an approximation.

The following example is covered by `Docs_Demo_IntegerBoundaries()` in
`src/docs_demos/base/docs_demo_standard_functions.cpp`, runnable through
`MML_DocsApp --integer-boundaries`. It also exercises Dirac approximations at
maximum positive `int` sharpness:

```cpp
const long long factorial = Functions::FactorialInt(20);
std::cout << factorial << '\n';  // 2432902008176640000
bool rejected = false;
try {
    Functions::FactorialInt(21);
}
catch (const DomainError&) {
    rejected = true;
    std::cout << "21! exceeds the exact signed 64-bit range\n";
}
const int n = std::numeric_limits<int>::max();
const DiracStep step(n);
const DiracSqr cauchy(n);
const Real x = Real(1) / static_cast<Real>(n);
std::cout << step(0) << '\n';
std::cout << cauchy(x) << '\n';  // n / (2*pi)
```

---

## Usage Patterns

### Direct Use
```cpp
using namespace MML::Functions;

Real result = Sin(Constants::PI / 2);
```

### In Function Objects
```cpp
#include <mml/base/Function.h>

// Create function object for sin(x)
auto f = RealFunction([](Real x) { 
    return Functions::Sin(x); 
});

Real val = f(Constants::PI / 2);  // 1.0
```

### In Algorithms
```cpp
// Numerical integration of special functions
Real integral = NumericalIntegration(
    [](Real x) { return Functions::Hermite(3, x) * Functions::Exp(-x*x); },
    -5.0, 5.0
);
```

---

## Notes

1. **Precision**: All functions operate on `Real` type (typically `double`)
2. **Domain**: Standard library domain restrictions apply (e.g., Asin requires |x| ≤ 1)
3. **Performance**: Inline wrappers have zero overhead over direct std:: calls
4. **Consistency**: Capitalized names distinguish from std:: versions

---

## See Also
- [Functions.md](Functions.md) - Function objects and composition
- [Interpolated_functions.md](Interpolated_functions.md) - Tabulated special functions
- [Polynoms.md](Polynoms.md) - Polynomial evaluation

---

## Runnable Examples

| Example | Description | Source |
|---------|-------------|--------|
| Standard Functions | Trig, special functions usage | See [`docs_demo_functions.cpp`](../../src/docs_demos/docs_demo_functions.cpp) |

**To run:** Build and execute `MML_DocsApp` target.
