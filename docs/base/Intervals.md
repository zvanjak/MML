# Intervals: Sets, Domains and Sampling

**Implementation:** [Intervals.h](../../mml/base/Intervals.h)

**Interface:** [IInterval.h](../../mml/interfaces/IInterval.h)

MML intervals represent sets of real values: finite ranges, rays, the real line,
periodically punctured ranges, and unions of pieces. This is **set algebra**, not
outward-rounded interval arithmetic for error enclosures. In particular, the
binary `-` operator means set difference, not subtraction of interval bounds.

All classes are **non-templated** and use the configured `MML::Real`.
`Interval` is an alias for `CompositeInterval`, not a separate implementation.
Use `ClosedInterval(0, 10)`, not a template specialization. Include the
implementation header and use the `MML` namespace; the examples also use
`<algorithm>`, `<cmath>`, `<iostream>`, `<limits>` and `<vector>`.

## Runnable Examples

Every C++ block below has a runnable counterpart in
[docs_demo_intervals.cpp](../../src/docs_demos/base/docs_demo_intervals.cpp).
The named example functions are called by `Docs_Demo_Intervals()`;
`Docs_Demo_IntervalContracts()` additionally checks endpoint predicates, hole
rejection, distinct budgets and finite windows. Failures throw rather than just
printing an expected result.

From the repository root on Windows:

```powershell
cmake --build build --config Release --target MML_DocsApp --parallel
& .\build\src\docs_demos\Release\MML_DocsApp.exe --interval-contracts
```

For a single-configuration Linux/macOS Release build, configure with
`cmake -B build -DCMAKE_BUILD_TYPE=Release`, build the same target, and run
`./build/src/docs_demos/MML_DocsApp --interval-contracts`.

## Contents

- [Types and endpoint semantics](#types-and-endpoint-semantics)
- [Set operations and predicates](#set-operations-and-predicates)
- [Compound properties and normalization](#compound-properties-and-normalization)
- [Function domains](#function-domains)
- [Chaining operations and using rays](#chaining-operations-and-using-rays)
- [Preparing root brackets](#preparing-root-brackets)
- [Sampling contracts](#sampling-contracts)
- [Error handling](#error-handling)
- [API inventory](#api-inventory)
- [Cost and numerical limitations](#cost-and-numerical-limitations)
- [Validation](#validation)

## Types and Endpoint Semantics

`IInterval` provides point membership, bounds, length/measure, a continuity flag
and sampling. `BaseInterval` implements the shared queries and sampling;
concrete subclasses supply membership semantics. `CompositeInterval` derives
directly from `IInterval` and adds set algebra and interval predicates.

| Concrete type / constructor | Mathematical set |
|---|---|
| `ClosedInterval(a, b)` | `[a, b]` |
| `OpenInterval(a, b)` | `(a, b)` |
| `OpenClosedInterval(a, b)` | `(a, b]` |
| `ClosedOpenInterval(a, b)` | `[a, b)` |
| `ClosedToInfInterval(a)` | `[a, +infinity)` |
| `OpenToInfInterval(a)` | `(a, +infinity)` |
| `NegInfToClosedInterval(b)` | `(-infinity, b]` |
| `NegInfToOpenInterval(b)` | `(-infinity, b)` |
| `CompleteRInterval()` | The real line |
| `CompleteRWithReccuringPointHoles(hole0, period)` | Real line excluding `hole0 + n * period`, integer `n` |
| `ClosedIntervalWithReccuringPointHoles(a, b, hole0, period)` | `[a, b]` with the same periodic exclusions |

The public names spell **`Reccuring`** as shown. Ordinary finite endpoints are
compared exactly: no epsilon merges nearby bounds or closes an open endpoint.
A closed interval with equal bounds is a singleton; the other three bounded
endpoint combinations with equal bounds are empty.

```cpp
void IntervalTypesExample()
{
	const ClosedInterval closed(0, 10);
	const OpenInterval open(0, 10);
	const OpenClosedInterval openClosed(0, 10);
	const ClosedOpenInterval closedOpen(0, 10);
	const ClosedToInfInterval right(0);
	const OpenToInfInterval openRight(0);
	const NegInfToClosedInterval left(0);
	const NegInfToOpenInterval openLeft(0);
	const CompleteRInterval line;
	const ClosedInterval singleton(3, 3);
	const OpenInterval empty(3, 3);

	if (!closed.contains(0) || !closed.contains(10) ||
		open.contains(0) || open.contains(10) || !open.contains(5) ||
		openClosed.contains(0) || !openClosed.contains(10) ||
		!closedOpen.contains(0) || closedOpen.contains(10) ||
		!right.contains(0) || !right.contains(100) || openRight.contains(0) ||
		!openRight.contains(100) || !left.contains(0) || !left.contains(-100) ||
		openLeft.contains(0) || !openLeft.contains(-100) ||
		!line.contains(-100) || !line.contains(100) ||
		!singleton.contains(3) || empty.contains(3))
		throw InvalidStateError("Interval type example failed");
	std::cout << "Finite endpoints, rays, real line, singleton and empty interval verified.\n";
}
```

Use the ray classes instead of passing IEEE infinity as a finite bound.
Internally, unbounded endpoints have distinct endpoint tags and use
`+/-std::numeric_limits<Real>::max()` as scalar bound sentinels. Bounds alone do
not distinguish a finite extreme endpoint from an unbounded endpoint. Use
finite query points; IEEE infinities and NaNs are not a portable membership
protocol across the concrete classes.

## Set Operations and Predicates

For ordinary pieces, union, intersection, difference and complement return
normalized `Interval` values. A `BaseInterval` converts implicitly to a composite,
so a simple interval can be passed to the static and member operations.
For operators, put a composite on the left, as in `Interval(a) & b`.

| Operation | Static form | Member form | Operator |
|---|---|---|---|
| Union | `Interval::Union(a, b)` | `a.unionWith(b)` | `a \| b` |
| Intersection | `Interval::Intersection(a, b)` | `a.intersectWith(b)` | `a & b` |
| Difference | `Interval::Difference(a, b)` | `a.differenceWith(b)` | `a - b` |
| Complement in the real line | `Interval::Complement(a)` | `a.complement()` | `~a` |

`contains(Real)` is available through every interval type.
`contains(const BaseInterval&)` and `intersects(const BaseInterval&)` belong to
**`CompositeInterval`**, not `IInterval` or `BaseInterval`. Containment means the
argument minus the receiving set is empty; it includes self-containment, empty
subsets and joint coverage by several pieces. Intersection predicates test an
actual nonempty intersection, not touching hulls.

An intersection includes a selected boundary only if both operands include that
point. At coincident finite boundaries an open endpoint wins; an unrelated open
endpoint does not make every result boundary open. Difference and complement
invert inclusion at a cut: `[0,10] \ [3,7]` is `[0,3) union (7,10]`.

```cpp
void IntervalSetOperationsExample()
{
	const ClosedInterval sleeping(15, 20);
	const ClosedInterval working(18, 24);
	const Interval optimal = Interval::Intersection(sleeping, working);
	const Interval sleepOnly = Interval::Difference(sleeping, working);
	const Interval either = Interval::Union(sleeping, working);
	const Interval outside = Interval::Complement(sleeping);
	const Interval set(sleeping);

	if (!optimal.contains(ClosedInterval(18, 20)) || optimal.contains(21) ||
		!sleepOnly.contains(15) || sleepOnly.contains(18) ||
		!either.contains(ClosedInterval(15, 24)) ||
		!outside.contains(14) || outside.contains(15) || !outside.contains(21) ||
		!set.contains(sleeping) || !set.intersects(working) ||
		!(set & working).contains(18) || !(set - working).contains(15) ||
		!(set | working).contains(24) || !(~set).contains(14) ||
		!set.intersectWith(working).contains(18) ||
		!set.differenceWith(working).contains(15) ||
		!set.unionWith(working).contains(24) || !set.complement().contains(14))
		throw InvalidStateError("Interval set-operation example failed");

	const Interval touching = Interval::Intersection(OpenInterval(0, 1), OpenInterval(1, 2));
	std::vector<Real> points;
	touching.getEquidistantCovering(1, points);
	if (!points.empty() || Interval(OpenInterval(0, 1)).intersects(OpenInterval(1, 2)))
		throw InvalidStateError("Open intervals must not intersect at an excluded endpoint");
	std::cout << "Set operations: [18,20], [15,18), [15,24] and the complement of [15,20].\n";
}
```

There is no public `IsEmpty()` query. For two ordinary simple intervals,
`!Interval(a).intersects(b)` determines that their intersection is empty. The
example also samples the known-empty result and checks the empty vector.
In general an empty sample vector is not a mathematical emptiness certificate:
the default window may miss the set, or a narrow interval may contain no
representable samples.

### Periodic-hole policy

Periodic-hole membership excludes lattice coordinates within
`Precision::RoundoffZeroThreshold` of an integer. This numerical roundoff rule
is separate from exact ordinary endpoint comparisons.

Store a hole-bearing concrete piece with `AddInterval(concretePiece)` to retain
its membership and supported sampling. `Normalize()` is a no-op if any stored
piece is noncontinuous. Algebra, interval-containment predicates, intersection
predicates, and reconstruction through `BaseInterval` reject hole-bearing pieces
with `NotImplementedError`; they do not silently replace them with a hull.

The existing `Docs_Demo_IntervalContracts()` checks this rejection as well as
membership and sampling of `ClosedIntervalWithReccuringPointHoles(0, 4, 1, 2)`.

## Compound Properties and Normalization

`AddInterval` **copies the concrete piece** and returns a reference for chaining.
It does not normalize automatically. For ordinary pieces, `Normalize()` sorts,
removes empty degenerate pieces and merges overlaps or touching pieces when at
least one touching endpoint is closed. Two open endpoints leave the missing
point excluded.

- `getLength()` is the **hull span**, maximum stored upper bound minus minimum
  stored lower bound.
- `getMeasure()` is the **sum of stored piece lengths**. Overlaps are counted
  repeatedly until ordinary pieces are normalized. It is not a membership-based
  integration of periodic roundoff exclusion bands.
- An empty default composite reports zero bounds, length and measure. A singleton
  also has length zero, so length is not an emptiness test. Stored empty pieces
  can affect bounds until normalization.
- `CompositeInterval::isContinuous()` currently returns `false` even for one
  connected normalized piece. It is not a connectivity or emptiness test.
- Reset a composite by assignment from `Interval{}`. There are no public
  `Clear()`, `NumIntervals()` or piece-iteration methods.

```cpp
void IntervalCompoundPropertiesExample()
{
	Interval compound;
	compound.AddInterval(ClosedInterval(0, 5))
		.AddInterval(ClosedInterval(10, 15))
		.AddInterval(ClosedInterval(20, 25));
	if (compound.getLowerBound() != 0 || compound.getUpperBound() != 25 ||
		compound.getLength() != 25 || compound.getMeasure() != 15 ||
		!compound.contains(3) || compound.contains(7) || !compound.contains(12))
		throw InvalidStateError("Compound interval properties failed");

	Interval overlapping;
	overlapping.AddInterval(ClosedInterval(0, 2)).AddInterval(ClosedInterval(1, 3));
	if (overlapping.getMeasure() != 4)
		throw InvalidStateError("Stored-piece measure must count an overlap before normalization");
	overlapping.Normalize();
	if (overlapping.getMeasure() != 3 || overlapping.getLength() != 3 ||
		!overlapping.contains(ClosedInterval(0, 3)) || overlapping.isContinuous())
		throw InvalidStateError("Normalized composite properties failed");
	overlapping = Interval{};
	std::vector<Real> points;
	overlapping.getEquidistantCovering(1, points);
	if (!points.empty()) throw InvalidStateError("Reset interval is not empty");
	std::cout << "Disjoint measure=15, hull length=25; overlapping stored measure=4, normalized=3.\n";
}
```

## Function Domains

For tangent, the excluded points are `pi/2 + n*pi` for **all integer** `n`.
Use the periodic-hole class for that global domain. A finite union of open
branches describes only the branches added, not the whole tangent domain.
Membership does not evaluate the function or certify that every included input
will produce a well-conditioned numerical result.

```cpp
void IntervalDomainsExample()
{
	const Real pi = Constants::PI;
	const CompleteRWithReccuringPointHoles tangentDomain(pi / 2, pi);
	Interval threeBranches;
	threeBranches.AddInterval(OpenInterval(-pi / 2, pi / 2))
		.AddInterval(OpenInterval(pi / 2, 3 * pi / 2))
		.AddInterval(OpenInterval(3 * pi / 2, 5 * pi / 2));
	if (!tangentDomain.contains(0) || tangentDomain.contains(pi / 2) ||
		!tangentDomain.contains(-pi) || !threeBranches.contains(0) ||
		!threeBranches.contains(pi) || threeBranches.contains(pi / 2) ||
		threeBranches.contains(-pi))
		throw InvalidStateError("Tangent-domain example failed");
	std::cout << "Global periodic tangent domain and three finite open branches verified.\n";
}
```

## Chaining Operations and Using Rays

Composite results feed directly into subsequent operations. Infinite endpoint
tags are preserved; they are not finite bounds passed as IEEE infinity.

```cpp
void IntervalSetChainExample()
{
	const ClosedInterval universe(0, 100);
	const ClosedInterval exclude1(20, 30);
	const ClosedInterval exclude2(60, 70);
	const Interval step1 = Interval::Difference(universe, exclude1);
	const Interval remaining = Interval::Difference(step1, exclude2);
	if (!remaining.contains(0) || remaining.contains(20) || remaining.contains(30) ||
		!remaining.contains(40) || remaining.contains(60) || remaining.contains(70) ||
		!remaining.contains(100) || remaining.getMeasure() != 80)
		throw InvalidStateError("Chained interval difference failed");

	const ClosedToInfInterval ray(0);
	const Interval overlap = Interval::Intersection(ray, ClosedInterval(5, 10));
	const Interval cutRay = Interval::Difference(ray, ClosedInterval(5, 10));
	if (!overlap.contains(ClosedInterval(5, 10)) ||
		!cutRay.contains(0) || cutRay.contains(5) || cutRay.contains(10) || !cutRay.contains(11))
		throw InvalidStateError("Ray set operations failed");
	std::cout << "Chained difference: [0,20) U (30,60) U (70,100].\n";
}
```

## Preparing Root Brackets

The continuous polynomial `(x - 1) * (x - 3)` has roots at 1 and 3. This example
handles exact sampled zeros separately from strict sign changes, avoids a
potentially overflowing product of endpoint values, and does not subtract one
from an empty vector's size. Six samples hit both roots exactly; eight samples
produce two sign-change brackets. Each bracket is checked against the known roots.

This is **not a general root-isolation guarantee**: coarse samples can miss
even-multiplicity roots or several crossings inside one cell. A sign change only
implies a root when the function is continuous on that bracket. Sampling alone
does not establish continuity across a gap or singularity.

```cpp
void IntervalBracketingExample()
{
	const auto f = [](Real x) { return (x - 1) * (x - 3); };
	const ClosedInterval searchRange(0, 5);
	for (const int sampleCount : {6, 8}) {
		std::vector<Real> points;
		searchRange.getEquidistantCovering(sampleCount, points);
		std::vector<Real> exactRoots;
		std::vector<ClosedInterval> brackets;
		for (size_t i = 0; i < points.size(); ++i) {
			const Real value = f(points[i]);
			if (value == 0) exactRoots.push_back(points[i]);
			if (i == 0) continue;
			const Real previous = f(points[i - 1]);
			if (previous != 0 && value != 0 && std::signbit(previous) != std::signbit(value))
				brackets.emplace_back(points[i - 1], points[i]);
		}
		if (exactRoots.size() + brackets.size() != 2 ||
			(sampleCount == 6 && exactRoots != std::vector<Real>{1, 3}) ||
			(sampleCount == 8 && brackets.size() != 2))
			throw InvalidStateError("Root sampling missed the two roots in this example");
		for (const auto& bracket : brackets)
			if (bracket.contains(1) == bracket.contains(3))
				throw InvalidStateError("Expected exactly one known root in each bracket");
		std::cout << sampleCount << " samples: " << exactRoots.size()
			<< " exact roots, " << brackets.size() << " sign-change brackets.\n";
	}
}
```

## Sampling Contracts

Sampling through `IInterval`, `BaseInterval` or `CompositeInterval` returns
**at most N distinct, ascending finite members**:

- Nonpositive requests clear the output without inspecting unused window arguments.
- Empty clipped sets produce no points; a singleton produces at most one.
- A one-point request uses an overflow-safe midpoint where it is a member.
- Closed endpoints are included where the budget permits; open endpoints move
  inward to the next representable value. Narrow ranges may have fewer than N
  representable members.
- Ordinary intervals use an equidistant grid. Hole-bearing sets use checked
  alternatives inside grid cells, so their samples need not be equidistant.
- Composite allocation uses normalized clipped occupied spans, not hull gaps.
  A small budget may omit components. Singletons receive at most one point.
  `Docs_Demo_IntervalContracts()` checks that 100 samples over `[0,1] union
  [100,101]` give 50 per piece.
- Bounded sets use their own bounds. If any component is unbounded, the default
  window is `[-Precision::UnboundedIntervalSamplingLimit,
  Precision::UnboundedIntervalSamplingLimit]`, intersected with the actual set.
- The four-argument overload takes an explicit **finite closed window**:
  `getEquidistantCovering(N, points, lower, upper)`. Custom `IInterval`
  implementations need not support it; the default throws `NotImplementedError`
  for positive requests.
- Sampling failures leave the old output unchanged. Membership callback
  exceptions propagate; unresolved periodic membership or candidate searches
  throw `NumericalMethodError` rather than publish partial success.

```cpp
void IntervalSamplingExample()
{
	std::vector<Real> points;
	ClosedInterval(0, 10).getEquidistantCovering(5, points);
	if (points != std::vector<Real>{0, Real(2.5), 5, Real(7.5), 10})
		throw InvalidStateError("Closed covering grid failed");
	OpenInterval(0, 10).getEquidistantCovering(5, points);
	if (points.size() != 5 || points.front() != std::nextafter(Real(0), Real(10)) ||
		points.back() != std::nextafter(Real(10), Real(0)))
		throw InvalidStateError("Open covering endpoints failed");
	ClosedInterval(3, 3).getEquidistantCovering(100, points);
	if (points != std::vector<Real>{3}) throw InvalidStateError("Singleton covering failed");

	const Real outside = 2 * Precision::UnboundedIntervalSamplingLimit;
	const ClosedToInfInterval ray(outside);
	const IInterval& view = ray;
	view.getEquidistantCovering(5, points);
	if (!points.empty()) throw InvalidStateError("Far ray must miss the default window");
	view.getEquidistantCovering(5, points, outside, outside + 4);
	if (points.size() != 5 || points.front() != outside || points.back() != outside + 4)
		throw InvalidStateError("Explicit window covering failed");
	std::cout << "Closed/open grids, distinct singleton and explicit ray window verified.\n";
}
```

Choose a sample count for the application and check the returned size rather than
assuming it equals the requested count. No fixed count guarantees integration
accuracy or makes a discontinuous composite a single valid integration bracket.

## Error Handling

| Condition | Exception / behavior |
|---|---|
| Nonfinite finite bounds, hole origin/period or used window bounds | `NumericInputError` |
| Reversed finite bounds or used window | `ArgumentError` |
| Zero periodic-hole period | `ArgumentError` |
| Finite negative period | Accepted; describes the same lattice as its positive counterpart |
| Unsupported hole-bearing algebra or abstract reconstruction | `NotImplementedError` |
| Positive window request on a custom implementation without that capability | `NotImplementedError` |
| Unrepresentable periodic lattice coordinate or unresolved nontrivial sampling cell | `NumericalMethodError` |
| Nonpositive sampling count | Clears output; no error, unused window arguments ignored |

Catch only errors the caller can handle. The deliberate failures below verify the
actual exception types. Use a ray constructor, not the rejected finite-bound form.

```cpp
void IntervalErrorsExample()
{
	bool reversedRejected = false;
	try {
		(void)ClosedInterval(10, 5);
	} catch (const ArgumentError& error) {
		reversedRejected = true;
		std::cout << error.what() << "\n";
	}
	bool nonfiniteRejected = false;
	try {
		(void)ClosedInterval(0, std::numeric_limits<Real>::infinity());
	} catch (const NumericInputError& error) {
		nonfiniteRejected = true;
		std::cout << error.what() << "\n";
	}
	if (!reversedRejected || !nonfiniteRejected)
		throw InvalidStateError("Invalid interval definitions were accepted");
}
```

## API Inventory

The tables below describe the current non-templated classes, not replacement
class declarations. Return types and `const` qualification matter: sampling
returns `void` and writes through the vector reference.

### `IInterval`

| Member | Contract |
|---|---|
| `virtual Real getLowerBound() const = 0` | Lower scalar bound |
| `virtual Real getUpperBound() const = 0` | Upper scalar bound |
| `virtual Real getLength() const = 0` | Hull span |
| `virtual Real getMeasure() const` | Defaults to `getLength()` |
| `virtual bool isContinuous() const = 0` | Implementation-specific flag; see composite caveat above |
| `virtual bool contains(Real x) const = 0` | Point membership |
| `virtual void getEquidistantCovering(int n, std::vector<Real>& points) const = 0` | Required sampling implementation |
| `virtual void getEquidistantCovering(int n, std::vector<Real>& points, Real lower, Real upper) const` | Optional explicit-window capability |
| `void GetEquidistantCovering(int n, std::vector<Real>& points) const` | Deprecated, nonvirtual two-argument forwarding alias; new code uses lowerCamel |
| `virtual ~IInterval()` | Polymorphic destruction |

Custom implementations override the lowerCamel name; the deprecated alias is
for callers, not an override point.

### `BaseInterval` and concrete pieces

`BaseInterval` has a protected constructor and remains abstract because point
membership is supplied by subclasses. It implements bounds, `getLength()`,
`isContinuous()` (true by default) and both sampling overloads; it inherits the
default `getMeasure()`. Periodic-hole subclasses return false for continuity.
The concrete public constructors are listed in the [types table](#types-and-endpoint-semantics).
There are no interval-containment or intersection predicates on `BaseInterval`.

### `CompositeInterval` / `Interval`

| Member | Result / behavior |
|---|---|
| `CompositeInterval()` | Empty stored piece list |
| `CompositeInterval(const BaseInterval& piece)` | Implicit reconstruction of an ordinary endpoint-typed piece; rejects holes |
| `template<class IntervalType> CompositeInterval& AddInterval(const IntervalType& piece)` | Copies a concrete `BaseInterval`-derived piece; chainable |
| `void Normalize()` | Canonicalizes ordinary pieces; no-op if any piece has holes |
| `CompositeInterval unionWith(const CompositeInterval& other) const` | Normalized union |
| `CompositeInterval intersectWith(const CompositeInterval& other) const` | Normalized intersection |
| `CompositeInterval differenceWith(const CompositeInterval& other) const` | Normalized difference |
| `CompositeInterval complement() const` | Complement in the real line |
| `bool contains(Real x) const` | Membership in any stored piece |
| `bool contains(const BaseInterval& other) const` | Whole simple interval is covered |
| `bool intersects(const BaseInterval& other) const` | Nonempty intersection with a simple interval |

The four static methods are `Union`, `Intersection` and `Difference`, each
returning `CompositeInterval` and accepting two `const CompositeInterval&`
arguments, and `Complement`, accepting one. The `|`, `&`, `-`, `~` operators
forward to the corresponding member methods. Ordinary simple arguments convert
implicitly; the predicate overloads themselves accept `BaseInterval`, not an
arbitrary composite.

The composite also overrides all required `IInterval` queries, `getMeasure()`
and both sampling overloads. These operations do not expose the stored vector.
There are no `RightRay`, `LeftRay`, `OpenRightRay`, `OpenLeftRay` or `RealLine`
aliases; use the concrete names in the types table.

## Cost and Numerical Limitations

These are source-level costs, not benchmark claims. With `n` and `m` stored
pieces and `k` generated candidate points:

- Composite point membership scans up to `n` pieces. Bounds and measure also
  scan stored pieces; `getLength()` performs two bound scans.
- Ordinary normalization sorts the stored pieces: `O(n log n)` time plus
  temporary storage.
- Union copies piece handles and normalizes; it is not an implicit side effect
  of `AddInterval`.
- Intersection checks all `n*m` pairs and normalizes the resulting pieces.
  Difference constructs a complement and then intersects. Containment uses
  difference, and the intersection predicate constructs an intersection result.
  These predicates are not promised constant-time or allocation-free.
- Sampling constructs and normalizes clipped allocation pieces, apportions
  budgets, checks candidate membership and sorts/deduplicates composite output.
  Hole searches add bounded refinement work dependent on `Real` precision.
  It is not accurately summarized as just a uniform-grid `O(k)` operation.
- For one ordinary piece per operand, set operations produce only a bounded
  number of pieces, but can still allocate.

Length/measure use ordinary `Real` subtraction and accumulation; extreme spans
can overflow, and unbounded sentinels are not mathematical measure values.
Sampling uses a separate overflow-safe finite-window path. Periodic-hole
membership can fail for an unrepresentable lattice coordinate. Neither these
queries nor sampling supply rigorous rounding-error enclosures.

## Validation

Executable examples live in
[docs_demo_intervals.cpp](../../src/docs_demos/base/docs_demo_intervals.cpp).
Interval tests are in [intervals_tests.cpp](../../tests/base/intervals_tests.cpp)
and [interval_set_contract_tests.cpp](../../tests/base/interval_set_contract_tests.cpp).
There are also degenerate-input and interface-default regressions. No fixed
assertion count or percentage coverage is claimed here.

For a focused Release double run on Windows, invoke the built executable
directly in one process:

```powershell
cmake --build build --config Release --target MML_Tests --parallel
& .\build\tests\Release\MML_Tests.exe '[intervals],[interval-contract]'
```

The repository's final code-change gate is one unfiltered direct Release
`MML_Tests` run; see [AGENTS.md](../../AGENTS.md). Builds enforce the configured
warning-as-error policy. Changing the example text requires updating and running
its demo counterpart, not merely making the Markdown look plausible.
