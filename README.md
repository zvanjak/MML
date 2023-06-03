# MML - Minimal Math Library
Basic mathematical functionality, contained in a single-header file.

**MML Vision**

For C++ developer
Who needs an easy to use math library to perform simple (and not so simple) math calculations
The Minimal Math Library
Is a single-header C++ library of classes and functions
That is trivial to use in any kind of project, is C++ 17 cross-platform compatible, and comes with a rich set of functionalities for working with vectors, matrices, linear systems, differential equations, coordinate systems and transformations.

**Math types**

- [vectors and matrices](/docs/Vectors_and_matrices.md)
  - Vector
  - VectorN<int N>
  - Matrix
  - MatrixNM<int N, int M>
- [tensors](/docs/Tensors.md) - TODO
  - Tensor2<int Dim>
  - Tensor3<int Dim>
  - Tensor4<int Dim>
- polynoms - TODO
  - Polynom
- [functions](/docs/Functions.md)
  - IRealFunction
  - IScalarFunction<int N>
  - IVectorFunction<int N>
  - InterpolatedFunction - TODO
- [geometry](/docs/Geometry.md) - TODO
  - Line2, Line3  
  - Plane
- [coordinate transformations](/docs/Coordinate_transformations.md)
  - spherical
  - cylindrical
  - general



**Algorithms**

- [linear alg. equations solvers](https://github.com/zvanjak/MML/blob/master/docs/Linear_equations_solvers.md)
  - GJ, LU, QR
- [numerical derivation](https://github.com/zvanjak/MML/blob/master/docs/Derivation.md)
- [numerical integration](https://github.com/zvanjak/MML/blob/master/docs/Integration.md)
- [vector operations](https://github.com/zvanjak/MML/blob/master/docs/Vector_operations.md)
  - grad, div
- [interpolation](https://github.com/zvanjak/MML/blob/master/docs/Interpolation.md)
- [differential equations solvers](https://github.com/zvanjak/MML/blob/master/docs/Differential_equations_solvers.md)
  - Runge-Kutta



