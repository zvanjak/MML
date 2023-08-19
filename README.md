# MML - Minimal Math Library
All your basic math needs, contained in a single-header file.

**MML Vision**
- For C++ developer, on Windows, Mac or Linux
- Who needs an easy to use math library to perform simple (and not so simple) math calculations
- The Minimal Math Library is a single-header C++ library of classes and functions
- That is trivial to use in any kind of project, is C++ 17 cross-platform compatible, and comes with a rich set of functionalities for working with vectors, matrices, linear systems, differential equations, coordinate systems and transformations.

**Intro**

Basic intro
- Based on Numerical recipes algorithms
- Pythonesque in its focus on simplicity of use (while trying as much as possible to retain C++ computational efficienvy)

**Basic math types**

- [Vectors](/docs/basic_types/Vector_types.md) - Vector, VectorN<int N>
- [Matrices](/docs/basic_types/Matrix_types.md) - Matrix, MatrixNM<int N, int M>
- [Tensors](/docs/basic_types/Tensors.md) - Tensor2<int Dim>, Tensor3<int Dim>, Tensor4<int Dim>
- [Polynoms](/docs/basic_types/Polynom.md) - Polynom, LinearFunctional, QuadraticForm
- [Functions](/docs/basic_types/Functions.md) - IRealFunction, IScalarFunction<int N>, IVectorFunction<int N>, InterpolatedFunction 1D, 2D, 3D
- [Geometry](/docs/basic_types/Geometry.md)
  - 2D geometry (Point2Cartesian, Point2Polar, Line2, SegmentLine2)
  - 3D geometry (Plane, Curves, Surfaces)
- [Coordinate transformations](/docs/basic_types/Coordinate_transformations.md)  
- [Algebra](/docs/Abstract_algebra.md) - groups, vector spaces, metric spaces
- [Differential geometry](/docs/Differential_geometry.md) - differential forms, manifolds

**Algorithms**

- [Linear alg. equations solvers](/docs/algorithms/Linear_equations_solvers.md) - GJ, LU, QR, SVD, Cholesky
- [Numerical derivation](/docs/algorithms/Derivation.md) - orders 1, 2, 4, 8
- [Numerical integration](/docs/algorithms/Integration.md) . Trapezoidal, Simpson, Romberg, Line & Surface integrals
- [Vector field operations](/docs/algorithms/Vector_field_operations.md) - grad, div, curl in cartesian, cylindrical and spherical coordinates
- [Interpolation](/docs/algorithms/Interpolation.md) - linear, polynomial, rational poly., spline (Bezier, B-spline, NURBS?)
- [Differential equations solvers](/docs/algorithms/Differential_equations_solvers.md) - Runge-Kutta
- [Root finding](/docs/algorithms/Root_finding.md)  
- [Function optimization](/docs/algorithms/Function_optimization.md)
- [Fourier transformation](/docs/algorithms/Fourier_transformation.md)
- [Statistics](/docs/algorithms/Statistics.md)

**Simple example**
~~~ c++
// vektori i matrice i funkcije? u interplayu 
// funkcija - deriviranje 
// lin system, zadan i riješen s dvi metode, usporedba rješenja
// ODE system
// interpolated as solution

    MML::Matrix<Real> m3(1, 3, {1.0, 1.0, 1.0});
    MML::Matrix<Real> m4(3, 4, {1.0, 0.0, 0.0, 0.0,
                                0.0, 1.0, 0.0, 0.0, 
                                0.0, 0.0, 1.0, 1.0});
    MML::Matrix<Real> m5 = m3 * m4;

    std::cout << "m3 = " << m3 << std::endl;
    std::cout << "m4 = " << m4 << std::endl;
    std::cout << "m3 * m4 = " << m5  << std::endl;
~~~

**Test beds**
- linear systems of equations
- real functions
- scalar functions
- vector functions
- ODE systems

**Precision testing**

**Speed testing**

