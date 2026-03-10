# Mathematical Background for PDE Solvers

> **Understanding the mathematics behind numerical PDE methods**

This document provides the mathematical foundations for the PDE solver module, including derivations, convergence theory, and key formulas.

---

## Table of Contents

1. [Partial Differential Equations Overview](#partial-differential-equations-overview)
2. [Elliptic PDEs](#elliptic-pdes)
3. [Parabolic PDEs](#parabolic-pdes)
4. [Finite Difference Methods](#finite-difference-methods)
5. [Sparse Linear Systems](#sparse-linear-systems)
6. [Iterative Solvers](#iterative-solvers)
7. [Preconditioners](#preconditioners)
8. [Convergence and Stability](#convergence-and-stability)
9. [Error Analysis](#error-analysis)

---

## Partial Differential Equations Overview

### Classification

PDEs are classified by the nature of their principal (highest-order) terms:

| Type | Canonical Form | Physical Example | Characteristics |
|------|---------------|------------------|-----------------|
| **Elliptic** | -∇²u = f | Steady-state heat, electrostatics | No time; equilibrium |
| **Parabolic** | ∂u/∂t = ∇²u | Heat diffusion | Time evolution; smoothing |
| **Hyperbolic** | ∂²u/∂t² = c²∇²u | Wave propagation | Time evolution; transport |

### Second-Order Classification

For a general second-order linear PDE in 2D:

$$Au_{xx} + Bu_{xy} + Cu_{yy} + Du_x + Eu_y + Fu = G$$

The discriminant $B^2 - 4AC$ determines the type:
- **Elliptic:** $B^2 - 4AC < 0$
- **Parabolic:** $B^2 - 4AC = 0$
- **Hyperbolic:** $B^2 - 4AC > 0$

### Boundary Conditions

| Type | Mathematical Form | Physical Meaning |
|------|-------------------|------------------|
| **Dirichlet** | $u = g$ on $\partial\Omega$ | Prescribed value |
| **Neumann** | $\frac{\partial u}{\partial n} = g$ on $\partial\Omega$ | Prescribed flux |
| **Robin (Mixed)** | $\alpha u + \beta\frac{\partial u}{\partial n} = g$ | Convective transfer |

---

## Elliptic PDEs

### Poisson Equation

The Poisson equation is the canonical elliptic PDE:

$$-\nabla^2 u = f \quad \text{in } \Omega$$

where the Laplacian operator is:

$$\nabla^2 u = \frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2} + \frac{\partial^2 u}{\partial z^2}$$

### Well-Posedness

For the Dirichlet problem $-\nabla^2 u = f$ with $u = g$ on $\partial\Omega$:

**Existence and Uniqueness:** A unique solution exists in $H^1(\Omega)$ if:
- $\Omega$ is a bounded domain with Lipschitz boundary
- $f \in L^2(\Omega)$ (or weaker: $f \in H^{-1}(\Omega)$)
- $g \in H^{1/2}(\partial\Omega)$

**Maximum Principle:** If $-\nabla^2 u \geq 0$ in $\Omega$, then:

$$\max_{\bar{\Omega}} u = \max_{\partial\Omega} u$$

### Variational Formulation

The Poisson problem has an equivalent variational (weak) form:

Find $u \in H^1_0(\Omega)$ such that for all $v \in H^1_0(\Omega)$:

$$\int_\Omega \nabla u \cdot \nabla v \, dx = \int_\Omega f v \, dx$$

Or equivalently, minimize the energy functional:

$$J(u) = \frac{1}{2}\int_\Omega |\nabla u|^2 \, dx - \int_\Omega f u \, dx$$

---

## Parabolic PDEs

### Heat Equation

The heat (diffusion) equation is the canonical parabolic PDE:

$$\frac{\partial u}{\partial t} = \alpha \nabla^2 u \quad \text{in } \Omega \times (0,T)$$

where $\alpha > 0$ is the thermal diffusivity.

### Initial-Boundary Value Problem

**Problem statement:**
$$\begin{align}
u_t &= \alpha \nabla^2 u + f & &\text{in } \Omega \times (0,T) \\
u &= g & &\text{on } \partial\Omega \times (0,T) \\
u(\cdot, 0) &= u_0 & &\text{in } \Omega
\end{align}$$

**Well-posedness:** Exists unique solution if $u_0 \in L^2(\Omega)$.

### Fundamental Solution

In unbounded domain $\mathbb{R}^d$, the fundamental solution is:

$$G(x,t) = \frac{1}{(4\pi\alpha t)^{d/2}} \exp\left(-\frac{|x|^2}{4\alpha t}\right)$$

Key properties:
- Integral over $\mathbb{R}^d$ equals 1
- Maximum at origin, decays exponentially
- Smoothing: infinite differentiability for $t > 0$

---

## Finite Difference Methods

### Taylor Series Basis

Finite difference approximations are derived from Taylor series:

$$u(x+h) = u(x) + h u'(x) + \frac{h^2}{2}u''(x) + \frac{h^3}{6}u'''(x) + O(h^4)$$

$$u(x-h) = u(x) - h u'(x) + \frac{h^2}{2}u''(x) - \frac{h^3}{6}u'''(x) + O(h^4)$$

### First Derivative Approximations

| Name | Formula | Order |
|------|---------|-------|
| Forward | $\frac{u(x+h) - u(x)}{h}$ | $O(h)$ |
| Backward | $\frac{u(x) - u(x-h)}{h}$ | $O(h)$ |
| Central | $\frac{u(x+h) - u(x-h)}{2h}$ | $O(h^2)$ |

### Second Derivative Approximation

Adding the Taylor expansions for $u(x+h)$ and $u(x-h)$:

$$u''(x) = \frac{u(x+h) - 2u(x) + u(x-h)}{h^2} + O(h^2)$$

This is the standard **second-order central difference** for the Laplacian.

### 2D Laplacian (5-Point Stencil)

$$\nabla^2 u \approx \frac{u_{i+1,j} + u_{i-1,j} + u_{i,j+1} + u_{i,j-1} - 4u_{i,j}}{h^2}$$

Stencil representation:

```
        [ 1 ]
    [1] [-4] [1]   × (1/h²)
        [ 1 ]
```

### 3D Laplacian (7-Point Stencil)

$$\nabla^2 u \approx \frac{u_{i+1,j,k} + u_{i-1,j,k} + u_{i,j+1,k} + u_{i,j-1,k} + u_{i,j,k+1} + u_{i,j,k-1} - 6u_{i,j,k}}{h^2}$$

---

## Sparse Linear Systems

### Discrete Poisson System

Discretizing $-\nabla^2 u = f$ with mesh size $h$ leads to:

$$\mathbf{A}\mathbf{u} = \mathbf{f}$$

where $\mathbf{A}$ is the **discrete Laplacian matrix**.

### 1D Discrete Laplacian

For n interior points:

$$\mathbf{A}_{1D} = \frac{1}{h^2}\begin{bmatrix}
2 & -1 & & & \\
-1 & 2 & -1 & & \\
& -1 & 2 & -1 & \\
& & \ddots & \ddots & \ddots \\
& & & -1 & 2
\end{bmatrix}$$

**Eigenvalues:**
$$\lambda_k = \frac{4}{h^2}\sin^2\left(\frac{k\pi h}{2}\right), \quad k = 1, 2, \ldots, n$$

**Condition number:**
$$\kappa(\mathbf{A}_{1D}) = \frac{\lambda_{\max}}{\lambda_{\min}} = \frac{\sin^2(n\pi h/2)}{\sin^2(\pi h/2)} \approx \frac{4}{\pi^2 h^2} = O(n^2)$$

### 2D Discrete Laplacian

For an $n \times n$ interior grid:

$$\mathbf{A}_{2D} = \mathbf{I} \otimes \mathbf{A}_{1D} + \mathbf{A}_{1D} \otimes \mathbf{I}$$

where $\otimes$ denotes Kronecker product.

**Eigenvalues:**
$$\lambda_{j,k} = \frac{4}{h^2}\left[\sin^2\left(\frac{j\pi h}{2}\right) + \sin^2\left(\frac{k\pi h}{2}\right)\right]$$

**Condition number:** $\kappa(\mathbf{A}_{2D}) = O(n^2) = O(1/h^2)$

### Properties of Discrete Laplacian

The discrete Laplacian matrix is:

| Property | Meaning | Implication |
|----------|---------|-------------|
| **Symmetric** | $A = A^T$ | Can use CG |
| **Positive Definite** | $x^TAx > 0$ for $x \neq 0$ | Unique solution exists |
| **Sparse** | $O(d)$ entries per row | Efficient SpMV |
| **M-matrix** | Positive diagonal, non-positive off-diagonal | Maximum principle holds |

---

## Iterative Solvers

### General Framework

Iterative methods generate a sequence $\mathbf{x}^{(k)}$ converging to $\mathbf{x}^*$ solving $\mathbf{Ax} = \mathbf{b}$.

**Convergence rate:** $\|\mathbf{e}^{(k+1)}\| \leq \rho \|\mathbf{e}^{(k)}\|$ where $\rho < 1$.

### Conjugate Gradient (CG)

For SPD matrix $\mathbf{A}$, CG minimizes:

$$\phi(\mathbf{x}) = \frac{1}{2}\mathbf{x}^T\mathbf{A}\mathbf{x} - \mathbf{b}^T\mathbf{x}$$

**Algorithm:**
```
r₀ = b - Ax₀
p₀ = r₀
for k = 0, 1, 2, ...
    αₖ = (rₖᵀrₖ) / (pₖᵀApₖ)
    xₖ₊₁ = xₖ + αₖpₖ
    rₖ₊₁ = rₖ - αₖApₖ
    if ‖rₖ₊₁‖ < tol then stop
    βₖ = (rₖ₊₁ᵀrₖ₊₁) / (rₖᵀrₖ)
    pₖ₊₁ = rₖ₊₁ + βₖpₖ
```

**Convergence bound:**

$$\|\mathbf{e}^{(k)}\|_A \leq 2\left(\frac{\sqrt{\kappa} - 1}{\sqrt{\kappa} + 1}\right)^k \|\mathbf{e}^{(0)}\|_A$$

For condition number $\kappa$, requires $O(\sqrt{\kappa})$ iterations for factor-of-2 reduction.

### BiCGSTAB

For non-symmetric systems, BiCGSTAB applies bi-conjugate gradients with stabilization.

**Cost per iteration:** 2 SpMV, 6 dot products, 6 axpy operations

**Convergence:** More erratic than CG; may stagnate for difficult problems.

### GMRES

Generalized Minimal Residual minimizes $\|\mathbf{b} - \mathbf{Ax}\|$ over Krylov subspace.

**Arnoldi process:** Builds orthonormal basis $\{v_1, v_2, \ldots, v_m\}$ for $\mathcal{K}_m(\mathbf{A}, \mathbf{r}_0)$.

**Memory:** Stores all basis vectors → $O(mn)$ for $m$ iterations.

**Restart:** After $m$ iterations, restart with current solution to limit memory.

---

## Preconditioners

### Preconditioning Concept

Instead of solving $\mathbf{Ax} = \mathbf{b}$, solve:

$$\mathbf{M}^{-1}\mathbf{Ax} = \mathbf{M}^{-1}\mathbf{b}$$

where $\mathbf{M} \approx \mathbf{A}$ but $\mathbf{M}^{-1}$ is cheap to apply.

**Goal:** $\kappa(\mathbf{M}^{-1}\mathbf{A}) \ll \kappa(\mathbf{A})$

### Jacobi Preconditioner

$$\mathbf{M} = \text{diag}(\mathbf{A})$$

**Properties:**
- Setup: $O(n)$
- Apply: $O(n)$ (element-wise division)
- Effectiveness: Reduces $\kappa$ by constant factor

For the discrete Laplacian: $M = (4/h^2)I$

### SSOR Preconditioner

Symmetric Successive Over-Relaxation:

$$\mathbf{M} = \frac{1}{\omega(2-\omega)}(\mathbf{D} + \omega\mathbf{L})\mathbf{D}^{-1}(\mathbf{D} + \omega\mathbf{U})$$

where $\mathbf{A} = \mathbf{L} + \mathbf{D} + \mathbf{U}$.

**Optimal ω for 2D Laplacian:**
$$\omega_{opt} = \frac{2}{1 + \sin(\pi h)} \approx 2 - 2\pi h$$

**Properties:**
- Setup: $O(n)$
- Apply: $O(nnz)$ (forward + backward sweep)
- Effectiveness: Reduces $\kappa$ to $O(1/h) = O(\sqrt{n})$

### ILU(0) Preconditioner

Incomplete LU with zero fill-in:

Compute $\mathbf{L}$ and $\mathbf{U}$ such that $\mathbf{A} \approx \mathbf{LU}$ where non-zeros are only in positions where $\mathbf{A}$ has non-zeros.

**Properties:**
- Setup: $O(nnz)$
- Apply: $O(nnz)$ (triangular solves)
- Effectiveness: Often 5-10× iteration reduction

### Preconditioner Comparison

| Preconditioner | $\kappa$ Reduction | Iterations for Laplacian |
|----------------|-------------------|--------------------------|
| None | 1× | $O(n)$ |
| Jacobi | ~2× | $O(n)$ |
| SSOR(optimal) | $O(\sqrt{n})$ | $O(n^{1/4})$ |
| ILU(0) | Variable | $O(n^{1/3})$ to $O(n^{1/2})$ |

---

## Convergence and Stability

### Consistency

A finite difference scheme is **consistent** if the truncation error vanishes as $h \to 0$:

$$\tau_h = L_h u - f_h \to 0 \quad \text{as } h \to 0$$

where $L_h$ is the discrete operator and $u$ is the exact solution.

### Stability

A scheme is **stable** if solutions remain bounded:

$$\|\mathbf{u}_h\| \leq C\|\mathbf{f}_h\|$$

for some constant $C$ independent of $h$.

**For elliptic problems:** Stability follows from positivity of the discrete operator.

### Convergence (Lax Equivalence Theorem)

For a consistent scheme:

$$\text{Stability} \iff \text{Convergence}$$

The error $\|\mathbf{u}_h - \mathbf{u}\|$ vanishes as $h \to 0$.

### CFL Condition for Parabolic PDEs

For the forward Euler discretization of $u_t = \alpha u_{xx}$:

$$u_j^{n+1} = u_j^n + \frac{\alpha\Delta t}{h^2}(u_{j+1}^n - 2u_j^n + u_{j-1}^n)$$

**Stability requires:**
$$\frac{\alpha\Delta t}{h^2} \leq \frac{1}{2}$$

Or equivalently:
$$\Delta t \leq \frac{h^2}{2\alpha}$$

**Multi-dimensional extension:**
$$\Delta t \leq \frac{h^2}{2d\alpha}$$

where $d$ is the number of spatial dimensions.

### Von Neumann Stability Analysis

For linear constant-coefficient problems, substitute $u_j^n = \xi^n e^{i k j h}$:

**Amplification factor:** $\xi = g(k, \Delta t, h)$

**Stability condition:** $|\xi| \leq 1$ for all wavenumbers $k$.

**Forward Euler heat equation:**
$$\xi = 1 - 4\frac{\alpha\Delta t}{h^2}\sin^2\left(\frac{kh}{2}\right)$$

Stability requires $\xi \geq -1$, giving the CFL condition.

---

## Error Analysis

### Truncation Error

The local truncation error measures how well the scheme approximates the PDE at a single point:

$$\tau_{i,j} = L_h[u]_{i,j} - f_{i,j}$$

For the 5-point Laplacian:
$$\tau_{i,j} = -\frac{h^2}{12}(u_{xxxx} + u_{yyyy}) + O(h^4)$$

### Global Error

The global error $e_h = u_h - u$ satisfies:

$$L_h e_h = -\tau_h$$

For stable schemes:
$$\|e_h\| \leq C\|\tau_h\| = O(h^p)$$

where $p$ is the order of accuracy.

### Convergence Rates

| Method | Spatial Order | Temporal Order | Overall |
|--------|---------------|----------------|---------|
| 2nd-order FD | $O(h^2)$ | — | $O(h^2)$ |
| Forward Euler | $O(h^2)$ | $O(\Delta t)$ | $O(h^2 + \Delta t)$ |
| Backward Euler | $O(h^2)$ | $O(\Delta t)$ | $O(h^2 + \Delta t)$ |
| Crank-Nicolson | $O(h^2)$ | $O(\Delta t^2)$ | $O(h^2 + \Delta t^2)$ |

### Richardson Extrapolation

If the error has form $e_h = Ch^p + O(h^{p+1})$, compute solutions on grids $h$ and $2h$:

$$u^{(ext)} = \frac{2^p u_h - u_{2h}}{2^p - 1}$$

This gives $O(h^{p+1})$ accuracy.

### Error Norms

| Norm | Formula | Physical Meaning |
|------|---------|------------------|
| **L∞ (max)** | $\max_i |e_i|$ | Worst-case error |
| **L² (RMS)** | $\sqrt{h^d \sum_i |e_i|^2}$ | Average error |
| **L¹** | $h^d \sum_i |e_i|$ | Total absolute error |

### Grid Convergence Study

To verify implementation, compute error on successively refined grids:

| Grid | h | $\|e\|$ | Rate |
|------|---|---------|------|
| $10$ | 0.1 | 1.23e-3 | — |
| $20$ | 0.05 | 3.08e-4 | 2.0 |
| $40$ | 0.025 | 7.71e-5 | 2.0 |
| $80$ | 0.0125 | 1.93e-5 | 2.0 |

**Rate calculation:**
$$p = \frac{\log(e_{2h}/e_h)}{\log(2)}$$

Expected: $p \approx 2$ for second-order methods.

---

## Key Formulas Reference

### Discrete Operators

| Operator | 1D | 2D | 3D |
|----------|----|----|-----|
| **Laplacian** | $\frac{u_{i-1} - 2u_i + u_{i+1}}{h^2}$ | 5-point | 7-point |
| **Gradient** | $\frac{u_{i+1} - u_{i-1}}{2h}$ | $(D_x u, D_y u)$ | $(D_x u, D_y u, D_z u)$ |

### Stability Criteria

| Scheme | Condition | Practical Limit |
|--------|-----------|-----------------|
| Forward Euler (1D) | $\nu = \alpha\Delta t/h^2 \leq 0.5$ | $\Delta t \leq h^2/(2\alpha)$ |
| Forward Euler (2D) | $\nu \leq 0.25$ | $\Delta t \leq h^2/(4\alpha)$ |
| Forward Euler (3D) | $\nu \leq 1/6$ | $\Delta t \leq h^2/(6\alpha)$ |
| Implicit schemes | Unconditional | Choose $\Delta t$ for accuracy |

### Condition Numbers

| Matrix | $\kappa$ | CG Iterations |
|--------|----------|---------------|
| 1D Laplacian ($n$ pts) | $O(n^2)$ | $O(n)$ |
| 2D Laplacian ($n \times n$) | $O(n^2)$ | $O(n)$ |
| Preconditioned (Jacobi) | $O(n^2)$ | $O(n)$ |
| Preconditioned (SSOR-opt) | $O(n)$ | $O(\sqrt{n})$ |

---

## Further Reading

### Textbooks

1. **LeVeque, R.J.** (2007). *Finite Difference Methods for Ordinary and Partial Differential Equations*. SIAM.
2. **Trefethen, L.N. & Bau, D.** (1997). *Numerical Linear Algebra*. SIAM.
3. **Saad, Y.** (2003). *Iterative Methods for Sparse Linear Systems*. SIAM.
4. **Strang, G.** (2007). *Computational Science and Engineering*. Wellesley-Cambridge.

### MML Documentation

- [API_Reference.md](API_Reference.md) - Complete class documentation
- [Iterative_Solvers.md](Iterative_Solvers.md) - Solver implementation details
- [Troubleshooting.md](Troubleshooting.md) - Common issues and solutions
- [Performance_Guide.md](Performance_Guide.md) - Optimization strategies
