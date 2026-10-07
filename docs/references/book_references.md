# Book References

This document lists the foundational textbooks and reference works that informed the design and implementation of MinimalMathLibrary.

> 📖 Book covers provided by [Open Library](https://openlibrary.org/) where available.

---

## 🔢 Numerical Methods & Scientific Computing

---

### 1. Numerical Recipes: The Art of Scientific Computing (3rd Edition)

<table>
<tr>
<td width="120" valign="top">
<img src="https://covers.openlibrary.org/b/isbn/9780521880688-M.jpg" alt="Numerical Recipes" width="100"/>
</td>
<td valign="top">

**Authors:** William H. Press, Saul A. Teukolsky, William T. Vetterling, Brian P. Flannery  
**Publisher:** Cambridge University Press (2007)  
**ISBN:** 978-0521880688

> *The* definitive reference for numerical algorithms. MML's implementation of many algorithms—including LU decomposition, QR factorization, ODE solvers (RK4, Dormand-Prince, Bulirsch-Stoer), root finding (Brent's method), optimization (Nelder-Mead), FFT, and random number generation—draws heavily from this masterwork.

**Key chapters used:**
- Ch. 2: Solution of Linear Algebraic Equations
- Ch. 4: Integration of Functions
- Ch. 5: Evaluation of Functions (Chebyshev, rational approximation)
- Ch. 7: Random Numbers
- Ch. 9: Root Finding and Nonlinear Sets of Equations
- Ch. 10: Minimization or Maximization of Functions
- Ch. 12: Fast Fourier Transform
- Ch. 17: Ordinary Differential Equations

</td>
</tr>
</table>

---

### 2. Matrix Computations (4th Edition)

<table>
<tr>
<td width="120" valign="top">
<img src="https://covers.openlibrary.org/b/isbn/9781421407944-M.jpg" alt="Matrix Computations" width="100"/>
</td>
<td valign="top">

**Authors:** Gene H. Golub, Charles F. Van Loan  
**Publisher:** Johns Hopkins University Press (2013)  
**ISBN:** 978-1421407944

> The authoritative treatise on numerical linear algebra. Essential for understanding the numerical stability and implementation details of matrix decompositions (LU, QR, SVD, Cholesky), eigenvalue algorithms, and iterative methods.

**Key topics used:**
- LU and Cholesky factorization with pivoting
- QR decomposition via Householder reflections
- Singular Value Decomposition (SVD)
- Eigenvalue computation for symmetric matrices (Jacobi, QR iteration)
- Condition numbers and error analysis

</td>
</tr>
</table>

---

### 3. Mathematical Methods for Physicists (7th Edition)

<table>
<tr>
<td width="120" valign="top">
<img src="https://covers.openlibrary.org/b/isbn/9780123846549-M.jpg" alt="Mathematical Methods for Physicists" width="100"/>
</td>
<td valign="top">

**Authors:** George B. Arfken, Hans J. Weber, Frank E. Harris  
**Publisher:** Academic Press / Elsevier (2013)  
**ISBN:** 978-0123846549

> Comprehensive coverage of mathematical methods essential for physics and engineering. Invaluable for vector calculus, coordinate systems, special functions, differential equations, and tensor analysis.

**Key chapters used:**
- Ch. 1-3: Vector Analysis, Coordinate Systems, Tensor Analysis
- Ch. 5: Infinite Series
- Ch. 7: Ordinary Differential Equations
- Ch. 9: Sturm-Liouville Theory
- Ch. 11-14: Special Functions (Bessel, Legendre, Hermite, Laguerre)
- Ch. 17: Integral Transforms (Fourier)

</td>
</tr>
</table>

---

### 4. Essential Mathematical Methods for the Physical Sciences

<table>
<tr>
<td width="120" valign="top">
<img src="https://covers.openlibrary.org/b/isbn/9780521761147-M.jpg" alt="Essential Mathematical Methods" width="100"/>
</td>
<td valign="top">

**Authors:** K. F. Riley, M. P. Hobson  
**Publisher:** Cambridge University Press (2011)  
**ISBN:** 978-0521761147

> Excellent companion to Arfken for applied mathematical methods. Clear treatment of vector calculus, differential equations, and integral theorems (Gauss, Stokes, Green).

**Key topics used:**
- Line, surface, and volume integrals
- Divergence theorem and Stokes' theorem verification
- Coordinate transformations
- Partial differential equations

</td>
</tr>
</table>

---

### 5. Numerical Methods for Engineers (7th Edition)

<table>
<tr>
<td width="120" valign="top">
<img src="https://covers.openlibrary.org/b/isbn/9780073397924-M.jpg" alt="Numerical Methods for Engineers" width="100"/>
</td>
<td valign="top">

**Authors:** Steven C. Chapra, Raymond P. Canale  
**Publisher:** McGraw-Hill (2015)  
**ISBN:** 978-0073397924

> Practical engineering-focused treatment of numerical methods. Excellent for understanding error analysis, interpolation, and the practical application of numerical algorithms.

**Key topics used:**
- Error propagation and numerical precision
- Interpolation methods (Newton, Lagrange, spline)
- Numerical differentiation and Richardson extrapolation
- Gauss quadrature derivation

</td>
</tr>
</table>

---

### 6. An Introduction to Numerical Analysis (2nd Edition)

<table>
<tr>
<td width="120" valign="top">
<img src="https://covers.openlibrary.org/b/isbn/9780471624899-M.jpg" alt="Introduction to Numerical Analysis" width="100"/>
</td>
<td valign="top">

**Author:** Kendall E. Atkinson  
**Publisher:** Wiley (1989)  
**ISBN:** 978-0471624899

> Rigorous mathematical treatment of numerical analysis with proofs and error bounds. Essential for understanding convergence properties and stability.

**Key topics used:**
- Polynomial interpolation error analysis
- Numerical integration error bounds
- Iterative methods for linear systems
- Nonlinear equation solvers convergence

</td>
</tr>
</table>

---

### 7. Applied Numerical Linear Algebra

<table>
<tr>
<td width="120" valign="top">
<img src="https://covers.openlibrary.org/b/isbn/9780898713893-M.jpg" alt="Applied Numerical Linear Algebra" width="100"/>
</td>
<td valign="top">

**Author:** James W. Demmel  
**Publisher:** SIAM (1997)  
**ISBN:** 978-0898713893

> Deep dive into numerical linear algebra with emphasis on algorithm design and floating-point considerations. Complements Golub & Van Loan.

**Key topics used:**
- Floating-point arithmetic and IEEE 754
- Perturbation theory and condition numbers
- Eigenvalue sensitivity analysis

</td>
</tr>
</table>

---

## 🌊 Differential Equations & Dynamical Systems

---

### 8. Solving Ordinary Differential Equations I: Nonstiff Problems (2nd Edition)

<table>
<tr>
<td width="120" valign="top">
<img src="https://covers.openlibrary.org/b/isbn/9783540566700-M.jpg" alt="Solving ODEs I" width="100"/>
</td>
<td valign="top">

**Authors:** Ernst Hairer, Syvert P. Nørsett, Gerhard Wanner  
**Publisher:** Springer (1993)  
**ISBN:** 978-3540566700

> The definitive reference for ODE solvers. MML's adaptive Runge-Kutta methods (Dormand-Prince 5(4), 8(5,3)) and embedded error estimators are based on this work.

**Key algorithms used:**
- Dormand-Prince coefficients
- Dense output formulas
- Step size control strategies
- Error estimation techniques

</td>
</tr>
</table>

---

### 9. Solving Ordinary Differential Equations II: Stiff and Differential-Algebraic Problems (2nd Edition)

<table>
<tr>
<td width="120" valign="top">
<img src="https://covers.openlibrary.org/b/isbn/9783540604525-M.jpg" alt="Solving ODEs II" width="100"/>
</td>
<td valign="top">

**Authors:** Ernst Hairer, Gerhard Wanner  
**Publisher:** Springer (1996)  
**ISBN:** 978-3540604525

> Essential for stiff equation solvers and implicit methods. Reference for Bulirsch-Stoer and implicit Runge-Kutta methods.

**Key topics used:**
- Stiffness detection
- Implicit methods and stability regions
- Bulirsch-Stoer algorithm with rational extrapolation

</td>
</tr>
</table>

---

### 10. Nonlinear Dynamics and Chaos (2nd Edition)

<table>
<tr>
<td width="120" valign="top">
<img src="https://covers.openlibrary.org/b/isbn/9780813349107-M.jpg" alt="Nonlinear Dynamics and Chaos" width="100"/>
</td>
<td valign="top">

**Author:** Steven H. Strogatz  
**Publisher:** Westview Press (2015)  
**ISBN:** 978-0813349107

> Beautifully written introduction to dynamical systems and chaos theory. Essential for the `Systems::` namespace—fixed point analysis, bifurcations, Lyapunov exponents.

**Key topics used:**
- Fixed point classification (nodes, spirals, saddles)
- Bifurcation diagrams (saddle-node, pitchfork, Hopf)
- Lyapunov exponents and chaos detection
- Lorenz, Rössler, and other classic systems

</td>
</tr>
</table>

---

## 📐 Differential Geometry & Vector Calculus

---

### 11. Differential Geometry of Curves and Surfaces

<table>
<tr>
<td width="120" valign="top">
<img src="https://covers.openlibrary.org/b/isbn/9780132125895-M.jpg" alt="Differential Geometry" width="100"/>
</td>
<td valign="top">

**Author:** Manfredo P. do Carmo  
**Publisher:** Prentice Hall (1976)  
**ISBN:** 978-0132125895

> Classic text on differential geometry. Foundation for parametric curve analysis—curvature, torsion, Frenet-Serret frames.

**Key topics used:**
- Curvature and torsion formulas
- Frenet-Serret apparatus
- Surface curvature (Gaussian, mean)
- First and second fundamental forms

</td>
</tr>
</table>

---

### 12. Div, Grad, Curl, and All That: An Informal Text on Vector Calculus (4th Edition)

<table>
<tr>
<td width="120" valign="top">
<img src="https://covers.openlibrary.org/b/isbn/9780393925166-M.jpg" alt="Div Grad Curl" width="100"/>
</td>
<td valign="top">

**Author:** H. M. Schey  
**Publisher:** W. W. Norton (2005)  
**ISBN:** 978-0393925166

> Accessible introduction to vector calculus. Excellent for understanding gradient, divergence, curl, and the integral theorems.

**Key topics used:**
- Physical interpretation of field operations
- Derivation of divergence and Stokes' theorems
- Coordinate system transformations

</td>
</tr>
</table>

---

## 📊 Statistics & Probability

---

### 13. Probability and Statistics for Engineers and Scientists (9th Edition)

<table>
<tr>
<td width="120" valign="top">
<img src="https://covers.openlibrary.org/b/isbn/9780321629111-M.jpg" alt="Probability and Statistics" width="100"/>
</td>
<td valign="top">

**Authors:** Ronald E. Walpole, Raymond H. Myers, Sharon L. Myers, Keying Ye  
**Publisher:** Pearson (2016)  
**ISBN:** 978-0321629111

> Standard reference for probability distributions and statistical methods. Used for the Statistics namespace.

**Key topics used:**
- Probability distributions (Normal, Binomial, Poisson, etc.)
- Hypothesis testing
- Regression analysis

</td>
</tr>
</table>

---

## 🧮 Special Functions & Approximation Theory

---

### 14. Handbook of Mathematical Functions

<table>
<tr>
<td width="120" valign="top">
<img src="https://covers.openlibrary.org/b/isbn/9780486612720-M.jpg" alt="Handbook of Mathematical Functions" width="100"/>
</td>
<td valign="top">

**Editors:** Milton Abramowitz, Irene A. Stegun  
**Publisher:** Dover Publications (1965)  
**ISBN:** 978-0486612720

> The classic reference for special functions and mathematical constants. Essential for Bessel functions, orthogonal polynomials, and error functions.

**Key sections used:**
- Gamma and Beta functions
- Error function and complementary error function
- Bessel functions
- Orthogonal polynomials (Legendre, Hermite, Laguerre, Chebyshev)

</td>
</tr>
</table>

---

### 15. NIST Digital Library of Mathematical Functions (DLMF)

<table>
<tr>
<td width="120" valign="top">
<img src="https://dlmf.nist.gov/style/DLMF-2x.png" alt="NIST DLMF" width="100"/>
</td>
<td valign="top">

**Editors:** Frank W. J. Olver, Daniel W. Lozier, Ronald F. Boisvert, Charles W. Clark  
**Publisher:** Cambridge University Press / NIST (2010)  
**URL:** https://dlmf.nist.gov/

> Modern successor to Abramowitz & Stegun. Online reference for special functions with updated formulas and computational methods.

</td>
</tr>
</table>

---

## 🔧 Software Engineering & C++ Implementation

---

### 16. The C++ Programming Language (4th Edition)

<table>
<tr>
<td width="120" valign="top">
<img src="https://covers.openlibrary.org/b/isbn/9780321563842-M.jpg" alt="C++ Programming Language" width="100"/>
</td>
<td valign="top">

**Author:** Bjarne Stroustrup  
**Publisher:** Addison-Wesley (2013)  
**ISBN:** 978-0321563842

> The authoritative C++ reference by the language creator. Essential for modern C++ design patterns used in MML.

</td>
</tr>
</table>

---

### 17. Effective Modern C++: 42 Specific Ways to Improve Your Use of C++11 and C++14

<table>
<tr>
<td width="120" valign="top">
<img src="https://covers.openlibrary.org/b/isbn/9781491903995-M.jpg" alt="Effective Modern C++" width="100"/>
</td>
<td valign="top">

**Author:** Scott Meyers  
**Publisher:** O'Reilly Media (2014)  
**ISBN:** 978-1491903995

> Best practices for modern C++. Informed MML's use of move semantics, smart pointers, lambdas, and template metaprogramming.

</td>
</tr>
</table>

---

## 📚 Additional References

*To be expanded with specific references for:*
- Quaternion algebra and 3D rotations
- Convex hull and computational geometry algorithms
- Delaunay triangulation and Voronoi diagrams
- Optimization algorithms (BFGS, conjugate gradient)
- Monte Carlo integration methods
- Symbolic computation

---

*Last updated: January 2026*
