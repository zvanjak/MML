# Function Spaces

Overview of function space abstractions and when to use them.

**Sources**: 
- `mml/interfaces/IFunction.h` - Function interfaces
- `mml/base/Function.h` - Function wrappers  
- `mml/core/OrthogonalBasis.h` - Orthogonal basis framework

---

## Overview

MML provides two complementary function space frameworks:

1. **Function Interfaces** - Type hierarchy for scalar, vector, and tensor-valued functions
2. **Orthogonal Bases** - Polynomial and trigonometric bases for function decomposition

MML 2.0 adds a third layer under `mml/core/FunctionSpaces/`: finite
function-space coordinates and operator discretization. This layer connects the
evaluation-oriented interfaces above to linear algebra:

```text
function f(x) -> coefficient vector c -> expansion u_N(x)
operator L   -> matrix A             -> algebraic system A c = b
```

For the roadmap that extends these pieces toward infinite-dimensional function
spaces, finite-dimensional projections, differential operators as matrices,
ODE/BVP solvers, and sparse large-`N` methods, see
[../Function_Spaces_Improvement.md](../Function_Spaces_Improvement.md).

---

## Current FunctionSpaces API

Include the aggregate header:

```cpp
#include <mml/core/FunctionSpaces.h>
using namespace MML;
using namespace MML::FunctionSpaces;
```

### Spaces, Trial Spaces, And Expansions

The first implemented types are one-dimensional:

| Type | Role |
|------|------|
| `FunctionSpace1D` | Abstract domain, weight, and inner-product descriptor. |
| `L2IntervalSpace` | Unweighted `L^2(a,b)` interval. |
| `WeightedL2IntervalSpace` | Weighted interval with user-provided weight function. |
| `TrialSpace1D` | Abstract finite trial space: dimension, basis values, optional nodes, derivative hooks. |
| `FunctionExpansion1D` | Coefficient vector plus trial-space association; evaluates `u_N(x)` and derivatives when available. |

FunctionSpaces follows the existing `VectorSpaces` idea: coordinates are plain
`Vector<Real>` values, and meaning comes from the associated space/basis.

### Projection And Interpolation

`OrthogonalBasisTrialSpace1D` wraps an existing `OrthogonalBasis`, so modal
projection can reuse the existing orthogonal-basis coefficient machinery:

```cpp
LegendreBasis basis;
OrthogonalBasisTrialSpace1D V(basis, 4);

auto expansion = ProjectL2(f, V);
Real value = expansion.evaluate(0.25);
```

Nodal trial spaces can use interpolation:

```cpp
ChebyshevCollocationSpace1D V(-1.0, 1.0, 20);
auto expansion = Interpolate(f, V);
```

Both projection and interpolation also have diagnostic wrappers:

```cpp
auto result = TryInterpolate(f, V);
if (result.success()) {
	Real value = result.expansion->evaluate(0.25);
}
```

### Chebyshev Collocation

`ChebyshevCollocationSpace1D(a,b,N)` provides left-to-right
Chebyshev-Lobatto nodes, Chebyshev basis values, and first/second
differentiation matrices. It is the first concrete spectral/collocation trial
space for differential operators.

```cpp
ChebyshevCollocationSpace1D V(-1.0, 1.0, 40);
Matrix<Real> D1 = V.firstDerivativeMatrix();
Matrix<Real> D2 = V.secondDerivativeMatrix();
```

### Operators, Boundary Conditions, And Dense BVPs

Linear 1D differential operators are represented by coefficient functions per
derivative order:

```cpp
auto L = LinearDifferentialOperator1D::SecondOrder(
	[](Real) { return -1.0; }, // u'' coefficient
	[](Real) { return  0.0; }, // u' coefficient
	[](Real) { return  0.0; }  // u coefficient
);
```

Boundary conditions are explicit descriptors:

```cpp
BoundaryConditions1D bc{
	BoundaryCondition1D::Dirichlet(-1.0, 0.0),
	BoundaryCondition1D::Dirichlet( 1.0, 0.0)
};
```

The dense BVP path builds a collocation matrix, applies boundary rows, solves
with existing dense linear solvers, and lifts the coefficient vector back into a
`FunctionExpansion1D`:

```cpp
auto result = SolveDenseCollocationBVP(V, L, rhs, bc);
if (result.success()) {
	Real u = result.solution->evaluate(0.5);
}
```

### Matrix-Free Operators

For large discretizations, FunctionSpaces also defines a matrix-free
`LinearOperator` interface with `rows()`, `cols()`, and `apply(...)`.
`DirichletSecondDerivativeOperator1D` is the first concrete finite-difference
stencil prototype.

---

## MML 2.0 FunctionSpaces Ownership

The implementation home for the new finite function-space layer is
`mml/core/FunctionSpaces/`, aggregated by `mml/core/FunctionSpaces.h`.

FunctionSpaces should be a bridge layer, not a duplicate linear algebra or
geometry subsystem. It owns the vocabulary for representing functions by finite
coordinates and for turning operators into algebraic systems:

- `FunctionSpace1D` and later higher-dimensional descriptors define domains,
	weights, and inner products.
- `TrialSpace1D` and concrete trial spaces define finite bases, nodes, and
	derivative capabilities.
- `FunctionExpansion1D` stores coefficients in an existing `Vector<Real>` and
	evaluates the represented function.
- Operator assembly maps trial-space coordinates into existing `Matrix<Real>`
	or sparse matrix types.

### Dependencies To Reuse

| Existing area | FunctionSpaces role |
|---------------|---------------------|
| `mml/core/VectorSpaces/` | Source of finite-coordinate vocabulary: bases, coordinate changes, linear maps, subspaces. |
| `mml/core/OrthogonalBasis.h` and `mml/core/OrthogonalBasis/` | Source of modal basis evaluation, weights, normalization, and coefficient projection ideas. |
| `mml/base/ChebyshevApproximation.h` | Existing concrete model for function -> coefficients -> evaluation. |
| `mml/base/Vector/` and `mml/base/Matrix/` | Storage for coefficients, nodal values, dense collocation and Galerkin matrices. |
| `mml/core/Integration/` | Inner products, projection integrals, and weak-form assembly. |
| `mml/base/SparseMatrix/` and `mml/core/SparseSolvers/` | Large-`N` sparse assembly and iterative solve path. |
| Typed points, tensors, metrics, and differential forms | Future tensor-product spaces, coordinate-aware fields, geometric PDEs, and weak-form extensions. |

### Non-Goals For The Initial Layer

- Do not create top-level `mml/functions/` or `mml/operators/` directories for
	the first implementation wave.
- Do not reimplement vector-space basis or linear-map machinery inside
	FunctionSpaces.
- Do not reimplement sparse matrices or iterative solvers; add only the glue or
	builder APIs that FunctionSpaces genuinely needs.
- Do not start with general finite elements or multi-dimensional PDE meshes;
	begin with 1D spectral/collocation and dense BVP workflows.

### Sparse Large-N Inventory

The sparse FunctionSpaces path should reuse the sparse infrastructure that is
already in core instead of introducing a parallel sparse layer.

| Existing type/header | Reuse for FunctionSpaces |
|----------------------|--------------------------|
| `mml/base/SparseMatrix/SparseMatrixCOO.h` | Incremental operator assembly through triplets, duplicate consolidation, row insertion, and existing stencil helpers. |
| `mml/base/SparseMatrix/SparseMatrixCSR.h` | Primary compute format for sparse matrix-vector products and iterative solvers. |
| `mml/base/SparseMatrix/SparseMatrixCSC.h` | Optional column-oriented format if later assembly or direct-solver workflows need it. |
| `mml/core/SparseSolvers/ConjugateGradient.h` | First solver for symmetric positive definite weak-form and Poisson-like systems. |
| `mml/core/SparseSolvers/BiCGSTAB.h` and `GMRES.h` | Nonsymmetric collocation, convection-diffusion, and general operator systems. |
| `mml/core/SparseSolvers/Preconditioners.h` | Reuse Jacobi/SSOR/ILU-style preconditioners before adding FunctionSpaces-specific policies. |

FunctionSpaces still needs only glue-level APIs:

- sparse operator assembly helpers that map local stencil, collocation, or weak
	form contributions into `SparseMatrixCOO<Real>`;
- conversion points from assembled COO to `SparseMatrixCSR<Real>`;
- adapter helpers between `Vector<Real>` coefficients and the sparse solver
	`std::vector<Real>` interface where needed;
- result mapping from sparse solver statuses into FunctionSpaces BVP/operator
	diagnostics.

The dense path should remain independent and usable without sparse headers.

### Sparse Assembly Interface Sketch

FunctionSpaces sparse assembly should use `SparseMatrixCOO<Real>` as the local
builder and `SparseMatrixCSR<Real>` as the compute handoff:

```cpp
SparseMatrixCOO<Real> coo(space.dimension(), space.dimension());

// Local row/stencil contribution.
coo.addEntry(row, column, value);
coo.addRow(row, columns, values);

// After all local contributions and boundary rows are inserted.
SparseMatrixCSR<Real> A(coo);
```

The first sparse FunctionSpaces API should be a small glue layer, for example:

```cpp
SparseMatrixCOO<Real> BuildSparseFiniteDifferenceOperator(
		const FiniteDifferenceGrid1D& grid,
		const LinearDifferentialOperator1D& op,
		const BoundaryConditions1D& bc);

SparseMatrixCSR<Real> ToCSR(SparseMatrixCOO<Real>& assembled);
```

Design constraints:

- local assembly writes triplets or whole rows into COO;
- duplicate contributions are intentionally allowed and consolidated by COO;
- boundary conditions should be inserted through explicit row replacement or
	row assembly helpers, mirroring the dense `ApplyBoundaryRow` path;
- solver-facing APIs should convert `Vector<Real>` to `std::vector<Real>` only
	at the sparse-solver boundary;
- solver statuses from `SparseSolvers::SolverResult<Real>` should map into
	FunctionSpaces BVP diagnostics rather than leaking package-specific messages.

### Sturm-Liouville And Modal Eigenproblem Plan

FunctionSpaces should treat Sturm-Liouville problems as generalized operator
eigenproblems:

```text
L[u] = lambda M[u]
A c = lambda B c
```

where `A` is the assembled stiffness/operator matrix, `B` is a mass or weight
matrix, and `c` are the coefficients of a `FunctionExpansion1D`.

Recommended staged API:

```cpp
struct GeneralizedEigenproblem1D {
		Matrix<Real> stiffness;
		Matrix<Real> mass;
		TrialSpace1D const* trial_space;
};

GeneralizedEigenproblem1D BuildSturmLiouvilleProblem(
		const TrialSpace1D& V,
		const IRealFunction& p,
		const IRealFunction& q,
		const IRealFunction& w,
		const BoundaryConditions1D& bc);
```

Validation examples:

- vibrating string: `-u'' = lambda u`, `u(0)=u(1)=0`, eigenvalues
	`lambda_n = n^2*pi^2`;
- weighted Sturm-Liouville problem with nontrivial `w(x)` once generalized mass
	matrices are stable;
- beam and Schrödinger examples after fourth-order and potential terms are
	supported.

Dependencies and gaps:

- dense generalized eigensolver support should be confirmed before exposing a
	high-level `SolveSturmLiouville` convenience API;
- sparse generalized eigensolvers are a later large-`N` feature;
- boundary-condition application should reuse the explicit dense/sparse row
	replacement helpers instead of embedding boundary hacks in eigen solvers;
- eigenvectors should lift back to `FunctionExpansion1D` values so mode shapes
	can be evaluated and visualized.

### Tensor-Product Spaces And Separable PDE Plan

Tensor-product FunctionSpaces are the path from 1D operators to separable 2D/3D
PDEs without immediately building a full finite-element mesh stack.

For two 1D trial spaces `Vx` and `Vy`, a tensor-product basis has functions

```text
phi_ij(x,y) = phi_i(x) psi_j(y)
```

and coefficients can use lexicographic indexing:

```cpp
int index(int i, int j) { return i + nx * j; }
```

Planned core types:

```cpp
class TensorProductTrialSpace2D;
class FunctionExpansion2D;

Matrix<Real> KroneckerProduct(const Matrix<Real>& A, const Matrix<Real>& B);
Matrix<Real> KroneckerSum(const Matrix<Real>& Ax, const Matrix<Real>& Ay);
```

For separable operators, the standard pattern is:

```text
Delta_2D = Dxx kron Iy + Ix kron Dyy
```

Dense prototypes should use existing `Matrix<Real>` first. Sparse versions
should assemble directly into `SparseMatrixCOO<Real>` to avoid forming dense
Kronecker products for large grids.

Relationship to existing tensor/forms work:

- `Tensor`, rank-1 typed vectors, covectors, metrics, and differential forms
	provide geometric meaning for coordinate-aware PDEs;
- FunctionSpaces should initially treat tensor-product spaces as algebraic
	products of 1D trial spaces;
- geometric PDEs can later attach frames, metrics, and forms without changing
	the coefficient-vector storage model.

Validation examples:

- separable Poisson problem on a rectangle with known sine-product solution;
- 2D Laplacian sparsity-pattern checks;
- equivalence between dense Kronecker assembly and sparse local assembly on
	small grids;
- mode-shape visualization from tensor-product eigenvectors.

---

## Function Interface Selection Guide

| Interface | Signature | Use Case |
|-----------|-----------|----------|
| `IRealFunction` | f: ℝ → ℝ | 1D integration, derivation, root finding |
| `IScalarFunction<N>` | f: ℝⁿ → ℝ | Scalar fields, gradient, multidimensional integration |
| `IVectorFunction<N>` | f: ℝⁿ → ℝⁿ | Vector fields, curl, divergence, Jacobians |
| `IVectorFunctionNM<N,M>` | f: ℝⁿ → ℝᵐ | General mappings, coordinate transforms |
| `IRealToVectorFunction<N>` | f: ℝ → ℝⁿ | Parametric curves (base class) |
| `IParametricCurve<N>` | r(t): ℝ → ℝⁿ | Curves with parameter bounds |
| `IParametricSurface<N>` | r(u,w): ℝ² → ℝⁿ | Parametric surfaces |
| `ITensorField2-5<N>` | f: ℝⁿ → Tensor | Higher-order tensor fields |

## Quick Reference

- **Real functions**: `IRealFunction` for scalar functions of one variable; pair with 1D integration/derivation.
- **Scalar fields**: `IScalarFunction<N>` for `f: ℝⁿ → ℝ`; use gradient/Laplacian (with metric for general coords).
- **Vector fields**: `IVectorFunction<N>` for `f: ℝⁿ → ℝⁿ`; use curl/divergence and Jacobians.
- **General mappings**: `IVectorFunctionNM<N,M>` for `f: ℝⁿ → ℝᵐ`; use with coordinate transformations.
- **Parametric curves/surfaces**: `IParametricCurve<N>`, `IParametricSurface<N>` for geometric modeling and ODE solutions.
- **Tensor fields**: `ITensorField<N>` for higher-order field modeling; combine with coordinate transforms and metric tensors.

---

## Orthogonal Basis Framework

MML provides a complete framework for orthogonal function bases with the abstract `OrthogonalBasis` class:

| Basis | Class | Domain | Weight w(x) | Key Application |
|-------|-------|--------|-------------|-----------------|
| Fourier | `FourierBasis` | [-L, L] | 1 | Periodic phenomena, signal processing |
| Legendre | `LegendreBasis` | [-1, 1] | 1 | Spherical problems, multipoles |
| Hermite | `HermiteBasis` | (-∞,∞) | e^(-x²) | Quantum oscillator, Gaussian processes |
| Chebyshev | `ChebyshevBasis` | [-1, 1] | 1/√(1-x²) | Polynomial approximation, spectral methods |
| Laguerre | `LaguerreBasis` | [0, ∞) | e^(-x) | Hydrogen atom, semi-infinite problems |

**→ See [OrthogonalBases.md](OrthogonalBases.md) for comprehensive documentation including:**
- Detailed API for each basis (Evaluate, WeightFunction, Normalization, Recurrence)
- Mathematical properties (recurrence relations, special values)
- Physical applications (quantum mechanics, signal processing, PDEs)
- Integration with Gaussian quadrature
- Code examples for function expansion and reconstruction

---

## Cross-Links

- [Functions.md](Functions.md) - Detailed function interface documentation
- [OrthogonalBases.md](OrthogonalBases.md) - **Comprehensive orthogonal basis documentation**
- [Vector_field_operations.md](Vector_field_operations.md) - Field differential operators
- [Curves_and_surfaces.md](Curves_and_surfaces.md) - Parametric geometry
- [Integration.md](Integration.md) - Gaussian quadrature using polynomial zeros