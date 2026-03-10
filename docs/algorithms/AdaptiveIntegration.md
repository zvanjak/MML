# Adaptive 2D/3D Integration

This document describes the adaptive integration algorithms for 2D and 3D integrals
using quadtree and octree subdivision with error-driven cell refinement.

## Overview

Adaptive integration methods subdivide the integration domain intelligently, concentrating
computational effort where the integrand is difficult (peaks, rapid changes) and using
fewer evaluations in smooth regions.

**Key Features:**
- Quadtree-based 2D integration with local error estimation
- Octree-based 3D integration with local error estimation  
- Tensor-product Gauss-Kronrod (GK15) error estimation
- Automatic cell subdivision where error exceeds tolerance
- Configurable tolerances, depth limits, and evaluation budgets

## API Reference

### 2D Integration

```cpp
#include "MML.h"
using namespace MML::Integration;

// Simple API
auto result = IntegrateAdaptive2D(
    [](Real x, Real y) { return x * y; },  // Function to integrate
    0.0, 1.0,                               // x bounds [x1, x2]
    0.0, 1.0,                               // y bounds [y1, y2]
    1e-8                                    // Tolerance (optional, default 1e-8)
);

// Config-based API for full control
AdaptiveConfig2D config;
config.absoluteTolerance(1e-10)
      .relativeTolerance(1e-8)
      .maxDepth(20)
      .maxEvaluations(1000000);

auto result = IntegrateAdaptive2D(f, x1, x2, y1, y2, config);
```

### 3D Integration

```cpp
auto result = IntegrateAdaptive3D(
    [](Real x, Real y, Real z) { return x * y * z; },
    0.0, 1.0,     // x bounds
    0.0, 1.0,     // y bounds
    0.0, 1.0,     // z bounds
    1e-8          // Tolerance
);
```

### Result Structure

Both `AdaptiveResult2D` and `AdaptiveResult3D` contain:

| Field | Type | Description |
|-------|------|-------------|
| `value` | Real | Integral estimate (sum over all cells) |
| `error_estimate` | Real | Accumulated error bound from all cells |
| `function_evaluations` | int | Total number of function evaluations |
| `cells_subdivided` | int | Number of quadtree/octree splits performed |
| `max_depth_reached` | int | Maximum subdivision depth reached |
| `converged` | bool | True if total error < requested tolerance |

## Algorithm Details

### Quadtree/Octree Subdivision

The algorithm proceeds as follows:

1. **Initial Cell**: Apply tensor-product GK15 integration to entire domain
2. **Error Check**: Compare cell error to tolerance × (cell_area / total_area)
3. **Accept or Subdivide**: 
   - If error ≤ scaled tolerance → accept cell result
   - Otherwise → subdivide into 4 (2D) or 8 (3D) children and recurse
4. **Accumulate**: Sum results from all converged leaf cells

```
    ┌─────────────────────────────┐
    │                             │
    │    Error too large?         │
    │    Subdivide into 4...      │
    │                             │
    ├──────────────┬──────────────┤
    │    ┌────┬───┐│              │
    │    │ ✓  │ ✓ ││              │
    │    ├────┼───┤│      ✓       │
    │    │ ✓  │ ✓ ││              │
    │    └────┴───┘│              │
    ├──────────────┼──────────────┤
    │              │              │
    │      ✓       │      ✓       │
    │              │              │
    └──────────────┴──────────────┘
```

### Tensor-Product GK15 Error Estimation

For 2D cells, we use full tensor-product Gauss-Kronrod quadrature:

- **15×15 = 225** evaluations per cell (GK15 on each axis)
- **Kronrod estimate**: Uses all 225 points with GK15 weights
- **Gauss estimate**: Uses embedded G7 nodes (7×7 = 49 points)
- **Error**: |Kronrod - Gauss| provides reliable error bound

For 3D cells:
- **15×15×15 = 3375** evaluations per cell
- **Error**: |K15³ - G7³| comparison

This approach provides proper multi-dimensional error estimation, unlike simpler
nested 1D integration which only gives per-dimension error bounds.

### Tolerance Scaling

Each cell's tolerance is scaled by its volume fraction:

```
cell_tolerance = max(abs_tol × (cell_volume / total_volume), 
                     rel_tol × |cell_value|)
```

This ensures:
- Small cells in subdivided regions have proportionally smaller tolerances
- The sum of all cell errors respects the total requested tolerance
- Relative tolerance handles functions with large values

## Examples

### Example 1: Basic Polynomial

```cpp
#include "MML.h"
using namespace MML::Integration;

int main() {
    // ∬ x*y dx dy over [0,1]² = 0.25
    auto f = [](Real x, Real y) { return x * y; };
    
    auto result = IntegrateAdaptive2D(f, 0.0, 1.0, 0.0, 1.0, 1e-10);
    
    std::cout << "Integral: " << result.value << "\n";           // 0.25
    std::cout << "Error: " << result.error_estimate << "\n";     // ~1e-15
    std::cout << "Evals: " << result.function_evaluations << "\n";
    std::cout << "Converged: " << result.converged << "\n";      // true
    
    return 0;
}
```

### Example 2: Peak Function (Adaptive Advantage)

```cpp
// Sharp Gaussian peak at (0.5, 0.5)
auto peak = [](Real x, Real y) {
    Real dx = x - 0.5;
    Real dy = y - 0.5;
    return std::exp(-100.0 * (dx*dx + dy*dy));
};

auto result = IntegrateAdaptive2D(peak, 0.0, 1.0, 0.0, 1.0, 1e-6);

// Adaptive concentrates evaluations near the peak
std::cout << "Cells subdivided: " << result.cells_subdivided << "\n";
std::cout << "Max depth: " << result.max_depth_reached << "\n";
```

For localized functions, adaptive integration typically uses **10-100x fewer** 
evaluations than uniform grids while achieving the same accuracy.

### Example 3: Divergence Theorem Verification

```cpp
// Verify ∭ ∇·F dV = ∬ F·n dS for F = (x, y, z)
// ∇·F = 3, so ∭ 3 dV over unit cube = 3

auto divergence = [](Real x, Real y, Real z) { return 3.0; };

auto result = IntegrateAdaptive3D(divergence, 
                                  0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1e-8);

std::cout << "Volume integral of divergence: " << result.value << "\n";  // 3.0
// This equals the surface flux (verified analytically)
```

### Example 4: Sphere Volume

```cpp
// Volume of unit sphere = 4π/3 ≈ 4.189
auto sphere_indicator = [](Real x, Real y, Real z) {
    return (x*x + y*y + z*z <= 1.0) ? 1.0 : 0.0;
};

auto result = IntegrateAdaptive3D(sphere_indicator,
                                  -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, 1e-3);

std::cout << "Sphere volume: " << result.value << "\n";  // ~4.189
```

Note: Discontinuous indicator functions are challenging and may require looser 
tolerances due to the discontinuity at the sphere boundary.

## When to Use Adaptive Integration

**Use adaptive integration when:**
- The integrand has localized features (peaks, near-singularities)
- Function evaluations are expensive and should be minimized
- You need reliable error estimates
- The domain is rectangular

**Consider alternatives when:**
- The integrand is very smooth (fixed-grid integration may suffice)
- Extreme precision is needed (may require custom subdivision strategies)
- The domain is non-rectangular (see variable-limit integration, future work)
- Very high dimensions (Monte Carlo methods scale better)

## Performance Characteristics

| Dimension | Evals/Cell | Smooth Function | Peaked Function |
|-----------|------------|-----------------|-----------------|
| 2D | 225 | 1-10 cells | 10-100+ cells |
| 3D | 3375 | 1-10 cells | 10-1000+ cells |

**Memory**: O(max_depth) stack depth for recursion

**Complexity**: O(n × cells) where n is evaluations per cell

## Known Limitations

1. **Rectangular domains only**: Current implementation requires rectangular 
   (2D) or box (3D) domains. Variable-limit integration is planned.

2. **Discontinuities**: Functions with sharp discontinuities may require 
   special handling or domain decomposition.

3. **3D evaluation cost**: Each 3D cell requires 3375 evaluations, which 
   can be expensive for deeply subdivided regions.

4. **Error underestimation**: While rare, the GK error estimate can 
   occasionally underestimate true error for pathological functions.

## References

1. Piessens, R., et al. (1983). *QUADPACK: A Subroutine Package for Automatic 
   Integration*. Springer-Verlag.

2. Press, W. H., et al. (2007). *Numerical Recipes: The Art of Scientific 
   Computing* (3rd ed.), Chapter 4.

3. Genz, A., & Malik, A. (1980). "Remarks on algorithm 006: An adaptive 
   algorithm for numerical integration over an N-dimensional rectangular 
   region". *Journal of Computational and Applied Mathematics*, 6(4), 295-302.

## See Also

- [GaussKronrod.h](../../../mml/core/Integration/GaussKronrod.h) - 1D Gauss-Kronrod implementation
- [Integration2DAdaptive.h](../../../mml/core/Integration/Integration2DAdaptive.h) - 2D implementation
- [Integration3DAdaptive.h](../../../mml/core/Integration/Integration3DAdaptive.h) - 3D implementation
