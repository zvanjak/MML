# Paper References

This document lists the academic papers and technical reports that informed specific algorithms and methods implemented in MinimalMathLibrary.

---

## 🌊 ODE Solvers

### 1. A Family of Embedded Runge-Kutta Formulae
**Authors:** J. R. Dormand, P. J. Prince  
**Journal:** Journal of Computational and Applied Mathematics, Vol. 6, No. 1, pp. 19-26 (1980)  
**DOI:** 10.1016/0771-050X(80)90013-3

> The foundational paper introducing the Dormand-Prince family of embedded Runge-Kutta methods. MML's `DormandPrince5_Stepper` (DOPRI5) implements the 5(4) method from this paper, using the 4th-order solution for error estimation.

**Algorithm implemented:** Dormand-Prince 5(4) with FSAL (First Same As Last) property

---

### Dense Output: Some Practical Runge-Kutta Formulas
**Author:** L. F. Shampine
**Journal:** Mathematics of Computation, Vol. 46, No. 173, pp. 135-150 (1986)
**DOI:** 10.1090/S0025-5718-1986-0815836-3

> Provides the quartic continuous extension used by MML's `DormandPrince5_Stepper`. Its coefficients combine the seven accepted-step stage derivatives without additional right-hand-side evaluations.

**Algorithm implemented:** Dormand-Prince 5(4) quartic dense output

---

### 2. High Order Embedded Runge-Kutta Formulae
**Authors:** P. J. Prince, J. R. Dormand  
**Journal:** Journal of Computational and Applied Mathematics, Vol. 7, No. 1, pp. 67-75 (1981)  
**DOI:** 10.1016/0771-050X(81)90010-3

> Historical background for high-order Dormand-Prince pairs. MML's `DormandPrince8_Stepper` preserves its public name but now implements Hairer's distinct DOP853 8(5,3) method.

**Related algorithm:** High-order Dormand-Prince embedded pairs

### DOP853 Implementation
**Authors:** E. Hairer, S. P. Norsett, G. Wanner
**Source:** Solving Ordinary Differential Equations I, Section II.5, and Hairer's original DOP853 software

> Defines the 12-stage 8(5,3) propagation method, blended error estimator, and seventh-order continuous extension implemented by `DormandPrince8_Stepper`.

**Algorithm implemented:** DOP853 8(5,3) with seventh-order dense output

---

### 3. A Variable Order Runge-Kutta Method for Initial Value Problems with Rapidly Varying Right-Hand Sides
**Authors:** J. R. Cash, A. H. Karp  
**Journal:** ACM Transactions on Mathematical Software, Vol. 16, No. 3, pp. 201-222 (1990)  
**DOI:** 10.1145/79505.79507

> Introduces the Cash-Karp variant of embedded Runge-Kutta with improved error estimation. Implemented as `CashKarp_Stepper` in MML.

**Algorithm implemented:** Cash-Karp 5(4) adaptive Runge-Kutta

---

### 4. Numerical Integration of Ordinary Differential Equations Based on Trigonometric Polynomials
**Author:** W. Gautschi  
**Journal:** Numerische Mathematik, Vol. 3, pp. 381-397 (1961)  
**DOI:** 10.1007/BF01386037

> Early work on exponential integrators for oscillatory problems. Background for understanding specialized ODE methods.

---

### 5. Solving Ordinary Differential Equations with Discontinuities
**Authors:** L. F. Shampine, S. Thompson  
**Journal:** Applied Numerical Mathematics, Vol. 28, pp. 49-63 (1998)

> Strategies for handling discontinuities and events in ODE integration. Informed MML's event detection capabilities.

---

## 🔢 Linear Algebra & Matrix Computation

### 6. An Analysis of the Total Least Squares Problem
**Authors:** G. H. Golub, C. F. Van Loan  
**Journal:** SIAM Journal on Numerical Analysis, Vol. 17, No. 6, pp. 883-893 (1980)  
**DOI:** 10.1137/0717073

> Foundation for total least squares and connections to SVD. Essential background for understanding `SVDecompositionSolver`.

---

### 7. Matrix Computations and the Condition Number
**Author:** J. H. Wilkinson  
**Journal:** The Computer Journal, Vol. 4, No. 4, pp. 281-286 (1962)

> Pioneering work on numerical stability and condition numbers. Fundamental for understanding why certain solvers are preferred for ill-conditioned systems.

---

### 8. The Rotation of Eigenvectors by a Perturbation
**Authors:** C. Davis, W. M. Kahan  
**Journal:** SIAM Journal on Numerical Analysis, Vol. 7, No. 1, pp. 1-46 (1970)  
**DOI:** 10.1137/0707001

> Analysis of eigenvalue sensitivity. Important for understanding the behavior of eigensolvers on nearly degenerate matrices.

---

## 🎯 Root Finding & Optimization

### 9. An Algorithm with Guaranteed Convergence for Finding a Zero of a Function
**Author:** R. P. Brent  
**Journal:** The Computer Journal, Vol. 14, No. 4, pp. 422-425 (1971)  
**DOI:** 10.1093/comjnl/14.4.422

> Introduction of Brent's method combining bisection, secant, and inverse quadratic interpolation. MML's `FindRootBrent` implements this algorithm.

**Algorithm implemented:** Brent's root-finding method

---

### 10. A Simplex Method for Function Minimization
**Authors:** J. A. Nelder, R. Mead  
**Journal:** The Computer Journal, Vol. 7, No. 4, pp. 308-313 (1965)  
**DOI:** 10.1093/comjnl/7.4.308

> The original Nelder-Mead simplex algorithm for derivative-free optimization. MML's `NelderMead<N>` optimizer implements this classic method.

**Algorithm implemented:** Nelder-Mead simplex optimization

---

### 11. A New Approach to Variable Metric Algorithms
**Author:** R. Fletcher  
**Journal:** The Computer Journal, Vol. 13, No. 3, pp. 317-322 (1970)  
**DOI:** 10.1093/comjnl/13.3.317

> Development of BFGS quasi-Newton method for optimization. Reference for advanced optimization methods.

---

## 📈 Numerical Integration

### 12. Gaussian Quadrature Formulas
**Author:** A. H. Stroud, D. Secrest  
**Publisher:** Prentice-Hall (1966)

> Comprehensive treatment of Gaussian quadrature with tables. Essential for implementing Gauss-Legendre, Gauss-Laguerre, and Gauss-Hermite integration.

---

### 13. A Note on the Relative Efficiency of Gauss-Kronrod Integration
**Author:** L. N. Trefethen  
**Journal:** BIT Numerical Mathematics, Vol. 45, No. 3, pp. 557-560 (2005)

> Analysis of Gauss-Kronrod efficiency. Background for MML's `IntegrateGaussKronrod` implementation.

---

### 14. Adaptive Numerical Integration by the Trapezoidal Rule
**Authors:** C. T. H. Baker, M. S. Sherwood  
**Journal:** Journal of the Institute of Mathematics and Its Applications, Vol. 19, pp. 67-78 (1977)

> Adaptive strategies for trapezoidal integration. Influenced MML's adaptive integration approach.

---

## 🦋 Dynamical Systems & Chaos

### 15. Determining Lyapunov Exponents from a Time Series
**Authors:** A. Wolf, J. B. Swift, H. L. Swinney, J. A. Vastano  
**Journal:** Physica D: Nonlinear Phenomena, Vol. 16, No. 3, pp. 285-317 (1985)  
**DOI:** 10.1016/0167-2789(85)90011-9

> The standard algorithm for computing Lyapunov exponents from time series and ODE systems. MML's `LyapunovAnalyzer::Compute` is based on this work.

**Algorithm implemented:** Wolf algorithm for Lyapunov spectrum computation

---

### 16. Ergodic Theory of Chaos and Strange Attractors
**Authors:** J.-P. Eckmann, D. Ruelle  
**Journal:** Reviews of Modern Physics, Vol. 57, No. 3, pp. 617-656 (1985)  
**DOI:** 10.1103/RevModPhys.57.617

> Foundational paper on chaos theory and attractors. Essential background for dynamical systems analysis.

---

### 17. Deterministic Nonperiodic Flow
**Author:** E. N. Lorenz  
**Journal:** Journal of the Atmospheric Sciences, Vol. 20, No. 2, pp. 130-141 (1963)  
**DOI:** 10.1175/1520-0469(1963)020<0130:DNF>2.0.CO;2

> The original paper introducing the Lorenz system. MML's `LorenzSystem` is based directly on this work.

**System implemented:** Lorenz equations (σ, ρ, β parameters)

---

### 18. An Equation for Continuous Chaos
**Author:** O. E. Rössler  
**Journal:** Physics Letters A, Vol. 57, No. 5, pp. 397-398 (1976)  
**DOI:** 10.1016/0375-9601(76)90101-8

> Introduction of the Rössler system. Implemented as `RosslerSystem` in MML.

**System implemented:** Rössler attractor

---

## 📐 Numerical Differentiation

### 19. Practical Extrapolation Methods: Theory and Applications
**Authors:** A. Sidi  
**Publisher:** Cambridge University Press (2003)  
**ISBN:** 978-0521661591

> Modern treatment of Richardson extrapolation. Informed MML's high-order derivative formulas (`NDer4`, `NDer6`, `NDer8`).

---

### 20. On the Computation of Derivatives of Functions of Several Variables
**Author:** A. Griewank  
**Journal:** SIAM Journal on Numerical Analysis, Vol. 26, No. 3, pp. 704-735 (1989)

> Foundational work on automatic differentiation and numerical gradient computation. Background for MML's gradient calculations.

---

## 📊 Signal Processing & FFT

### 21. An Algorithm for the Machine Calculation of Complex Fourier Series
**Authors:** J. W. Cooley, J. W. Tukey  
**Journal:** Mathematics of Computation, Vol. 19, No. 90, pp. 297-301 (1965)  
**DOI:** 10.2307/2003354

> The original FFT paper that revolutionized signal processing. MML's `FFT::Transform` implements variants of this algorithm.

**Algorithm implemented:** Cooley-Tukey FFT

---

### 22. Real-Valued Fast Fourier Transform Algorithms
**Authors:** H. V. Sorensen, D. L. Jones, M. T. Heideman, C. S. Burrus  
**Journal:** IEEE Transactions on Acoustics, Speech, and Signal Processing, Vol. 35, No. 6, pp. 849-863 (1987)

> Optimized algorithms for real-valued FFT. Informed MML's `FFT::RealFFT` implementation.

---

## 🔄 Coordinate Systems & Geometry

### 23. Quaternions and Rotation Sequences
**Author:** J. B. Kuipers  
**Publisher:** Princeton University Press (1999)  
**ISBN:** 978-0691102986

> Comprehensive treatment of quaternion algebra and 3D rotations. Foundation for MML's `Quaternion` class and rotation operations.

---

### 24. Geometric Tools for Computer Graphics
**Authors:** P. J. Schneider, D. H. Eberly  
**Publisher:** Morgan Kaufmann (2003)  
**ISBN:** 978-1558605947

> Practical algorithms for computational geometry. Reference for intersection tests, distance calculations, and geometric primitives in MML.

---

## 📚 Additional Papers to Include

*This list will be expanded with references for:*
- Delaunay triangulation algorithms
- Convex hull computation (Graham scan, Quickhull)
- KD-tree construction and nearest neighbor search
- Spline interpolation theory
- Monte Carlo integration error analysis
- Implicit ODE methods for stiff systems

---

*Last updated: January 2026*