# MML Precision Testing Documentation

## Overview

The MML Precision Testing Framework provides comprehensive verification of numerical algorithm accuracy across all major computational areas of the library. This documentation captures the methodology, findings, and insights from extensive precision testing.

## Quick Reference: Precision Summary

| Algorithm Category | Best Precision | Typical Precision | Notes |
|-------------------|----------------|-------------------|-------|
| **Derivation** | 10⁻¹⁴ | 10⁻¹⁰ to 10⁻¹² | Step size dependent |
| **Integration** | 10⁻¹⁵ | 10⁻¹⁰ to 10⁻¹⁴ | Method and subdivisions dependent |
| **Interpolation** | 10⁻¹⁶ | 10⁻¹² to 10⁻¹⁵ | Polynomial degree dependent |
| **ODE Solvers** | 10⁻¹⁴ | 10⁻⁸ to 10⁻¹² | Method order dependent |
| **Root Finding** | 10⁻¹⁵ | 10⁻¹⁴ to 10⁻¹⁵ | Brent/Ridders most reliable |
| **Linear Algebra** | 10⁻¹⁶ | 10⁻¹⁴ to 10⁻¹⁶ | Condition number dependent |

## Documentation Structure

```
docs/testing_precision/
├── README.md                           # This file - overview and navigation
├── FRAMEWORK.md                        # PrecisionTestFramework.h design and API
├── DERIVATION_ANALYSIS.md              # Numerical differentiation precision
├── INTEGRATION_ANALYSIS.md             # Numerical integration precision
├── INTERPOLATION_ANALYSIS.md           # Interpolation methods precision
├── ODE_SOLVER_ANALYSIS.md              # ODE solver precision and convergence
├── ROOT_FINDING_ANALYSIS.md            # Root finding algorithm comparison
├── LINEAR_ALGEBRA_ANALYSIS.md          # Matrix solver and decomposition precision
└── METHODOLOGY.md                      # Testing methodology and best practices
```

## Key Findings

### 1. Machine Precision Limits
- IEEE 754 double precision: ~15-16 significant decimal digits
- Achieved 10⁻¹⁶ accuracy in many well-conditioned problems
- Catastrophic cancellation limits accuracy in ill-conditioned problems

### 2. Algorithm Selection Guidelines

**For Derivatives:**
- Use **NDer6** or **NDer8** for highest accuracy (10⁻¹² to 10⁻¹⁴)
- **NDer4** provides good balance of speed and accuracy (10⁻⁸ to 10⁻¹⁰)
- Use **NDer2** (central differences) for simple cases
- Use **NDer1** (forward differences) only for boundary conditions

**For Integration:**
- Use **Romberg** for smooth functions (10⁻¹⁵ achievable)
- Use **Gauss-Legendre** for polynomial-like integrands
- Use **Simpson** for general purpose

**For Interpolation:**
- **Barycentric** interpolation for stability at high degrees
- **Cubic splines** for smooth interpolation with C² continuity
- Beware of Runge phenomenon at high polynomial degrees

**For ODEs:**
- Use **DP8** for high accuracy requirements
- Use **RK4** for good balance of speed and accuracy
- Use **Verlet/Leapfrog** for energy-conserving problems

**For Root Finding:**
- Use **Brent** for reliability (guaranteed convergence)
- Use **Newton** when derivative is available and fast convergence needed
- Use **Ridders** as a good general-purpose method

**For Linear Systems:**
- Use **Cholesky** for SPD matrices (most efficient)
- Use **QR** for better numerical stability
- Use **SVD** for ill-conditioned or rank-deficient problems

### 3. Common Pitfalls

1. **Step size selection** - Too small causes roundoff, too large causes truncation
2. **Ill-conditioning** - Check condition numbers before trusting results
3. **Convergence verification** - Always verify convergence for iterative methods
4. **Boundary effects** - Special handling needed at domain boundaries

## Running the Tests

```bash
# Build the precision testing application
cmake --build build --target MML_PrecisionTestingApp

# Run all tests
./build/src/testing_precision/Release/MML_PrecisionTestingApp

# Run specific category
./build/src/testing_precision/Release/MML_PrecisionTestingApp derivation
./build/src/testing_precision/Release/MML_PrecisionTestingApp integration
./build/src/testing_precision/Release/MML_PrecisionTestingApp interpolation
./build/src/testing_precision/Release/MML_PrecisionTestingApp ode
./build/src/testing_precision/Release/MML_PrecisionTestingApp roots
./build/src/testing_precision/Release/MML_PrecisionTestingApp linalg
```

## Framework Components

- **PrecisionTestFramework.h** - Core testing infrastructure
- **test_precision_derivation.cpp** - Derivative tests
- **test_precision_integration.cpp** - Integration tests
- **test_precision_interpolation.cpp** - Interpolation tests
- **test_precision_ode.cpp** - ODE solver tests
- **test_precision_roots.cpp** - Root finding tests
- **test_precision_linalg.cpp** - Linear algebra tests

## Related Documentation

- [Framework Design](FRAMEWORK.md) - Detailed API documentation
- [Methodology](METHODOLOGY.md) - How tests are designed and validated
- Individual analysis documents for each algorithm category

---

*Last updated: December 2024*
*MML Version: 1.0*
