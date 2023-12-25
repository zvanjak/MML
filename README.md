# MML - Minimal Math Library
All your basic (numerical) math needs, contained in a single-header file.

**MML Vision**
- For a C++ developer, on Windows, Mac or Linux
- Who needs an easy to use math library to perform simple (and not so simple) numerical math calculations
and is intimidated by the complexity of existing math libraries
- The Minimal Math Library is a single-header C++ library of classes and functions
general purpose, pythonesque, easy to use, and efficient
- That is trivial to use in any kind of project, is C++ 20 cross-platform compatible, and comes with a rich set of functionalities for working with vectors, matrices, tensors, polynoms, real, scalar and vector functions, linear systems, differential equations, coordinate systems and transformations, with algorithms for derivation, integration, root finding, interpolation, optimization, statistics, and more.

**Basic facts**
- As of now, and for forseable future, this is for personal and educational use only (see Licensing at the end)
- Single-header -> include and go 
- Cross-platform -> tested with Visual Studio Code on Windows (MSVC++, clang), Mac (gcc, clang), Linux (gcc)
- Pythonesque in its focus on simplicity of use, and focused on faitful modeling of mathematical entities (while trying as much as possible to retain C++ computational efficiency)

**Is there really a need for another C++ math library?**
- Main benefit, and the reason I did all this is attempt at completeness, and simplicity of use
- Not focused on efficiency so when going gets rough - Boost, Eigen, Armadillo, etc.

**Core**
- [Core utils](/docs/core/CoreUtils.md) - general utilities
- [Vectors](/docs/core/Vector_types.md) - Vector, VectorN<int N>
- [Matrices](/docs/core/Matrix_types.md) - Matrix, MatrixNM<int N, int M>, MatrixSymm, MatrixBandDiag, MatrixTridiag
- [Matrix utils](/docs/core/MatrixUtils.md) - general utilities
- [Algebra](/docs/core/Algebra.md) - groups
- [Vector spaces](/docs/core/Vector_spaces.md) - vector space, normed vector space, metric (Hilbert) space
- [Linear alg. equations solvers](/docs/core/Linear_equations_solvers.md) - GJ, LU, QR, SVD, Cholesky
- [Polynoms](/docs/core/Polynom.md) - Polynom
- [Function types](/docs/core/Function_types.md) - IRealFunction, IScalarFunction<int N>, IVectorFunction<int N>, IParametricCurve<N>, IParametricSurface<N>
- [Interpolated functions](/docs/core/Interpolated_functions.md) - linear, polynomial, rational poly., spline (Bezier, B-spline, NURBS?)
- [Field operations](/docs/core/Vector_field_operations.md) - grad, div, curl, Laplacian in general, Cartesian, cylindrical and spherical coordinates
- [Operators](/docs/core/Operators.md)  - linear functional, quadratic form, linear operator
- [Numerical derivation](/docs/core/Derivation.md) - orders 1, 2, 4, 6, 8
- [Numerical integration](/docs/core/Integration.md) - Trapezoidal, Simpson, Romberg algorithms
- [Geometry basic types](/docs/core/Geometry.md) - points, vectors in Cartesian, Polar, Cylindrical and Spherical coordinates
- [Tensors](/docs/core/Tensors.md)  - Tensor2<int Dim>, Tensor3<int Dim>, Tensor4<int Dim>
- [Metric tensor](/docs/core/Metric_tensor.md) - predefined metric tensors in Cartesian, Cylindrical and Spherical coordinates
- [Coordinate transformations](/docs/core/Coordinate_transformations.md) - General, Cartesian, Cylindrical, Spherical  

**Basic math types**
- [2D & 3D geometry](/docs/basic_types/Geometry_2D_3D.md) - Lines, Planes, Bodies
- [Curves & Surfaces](/docs/basic_types/Curves_and_surfaces.md) - predefined curves and surfaces
- [Dirac delta function](/docs/basic_types/Dirac_delta_function.md)- predefined distributions for representing Dirac delta function
- [Fields](/docs/basic_types/Fields.md) - predefined fields
- [Functions](/docs/basic_types/Functions.md) - definition of available standard functions
- [Function spaces](/docs/basic_types/Function_spaces.md) - Hermitian, Legendre, Laguerre, Chebyshev, Fourier spaces
- [ODE system](/docs/basic_types/ODE_system.md) - represents a system of ordinary differential equations

**Algorithms**
- [Eigen solvers](/docs/algorithms/Eigen_solvers.md) - solving eigenvalue problems for symmetric and non-symmetric matrices
- [Differential geometry](/docs/algorithms/Differential_geometry.md) - for curves only, so far
- [Function analyzer](/docs/algorithms/Function_analyzer.md) - analyzing functions, finding roots, extrema, inflection points, etc.
- [ODE system solvers](/docs/algorithms/Differential_equations_solvers.md) - solvers for systems of ordinary differential equations
- [Path integration](/docs/algorithms/Path_integration.md) - calculating line and work integrals
- [Multidim integration](/docs/algorithms/Multidim_integration.md) - calculating surface and volume integrals
- [Root finding](/docs/algorithms/Root_finding.md) - TODO
- [Statistics](/docs/algorithms/Statistics.md) - basics - avg, std, var, cov, corr
- [Fourier transformation](/docs/algorithms/Fourier_transformation.md) - TODO

**Test beds**
There is a multitude of test beds for each of the above algorithms, in the [testbeds](/docs/testbeds/testbeds_intro.md) folder.
They are used for verification of implemented algorithms, but can also be used as examples of how to use the library.

- real, scalar and vector functions
- linear systems of equations
- ODE systems
- parametric curves & surfaces

**Before intro, couple of real examples what it is for**
- primjer s transformacijama koordinata
- proracun sudara dva tijela
- proračun work integrala po zatvorenoj petlji u solenoidalnom polju?
- proračun momenta inercije za tijelo s diskretnim skupom masa

**Intro examples**

***Vectors, matrices***
~~~ c++
    // TODO - finish vectors, matrices example
    // matrica je operator, pa vektor se mnozi
    // dodati i rad real , complex matrice!
    Vector<double>  vec_dbl_3({ 1.0, 2.0, 3.0 }); 
    VectorComplex   vec_cmplx_2({ Complex(1,1), Complex(-1,2), Complex(2, -0.5) });

    Matrix<Real>   mat_3x3{ 3, 3, { 1.0, 2.0, -1.0, 
                                   -1.0, 5.0, 6.0, 
                                    3.0, 1.0, 1.0 }};  
    MatrixComplex  mat_cmplx(2,2, { Complex(1,1),  Complex(-1,2), 
                                    Complex(2, -0.5), Complex(1,1) });
    Matrix<Real>   unit_mat3 = MML::Matrix<Real>::GetUnitMatrix(3);
    // TODO - staviti link na "complete example with output"
~~~

***Solving linear systems of equations and calculating eigenvalues***
~~~ c++
 	Matrix<Real>    mat{5, 5, { 1.4, 2.1, 2.1, 7.4, 9.6,
                                1.6, 1.5, 1.1, 0.7, 5.0,
                                3.8, 8.0, 9.6, 5.4, 8.8,
                                4.6, 8.2, 8.4, 0.4, 8.0,
                                2.6, 2.9, 0.1, 9.6, 7.7 } };
	Vector<Real> 	rhs{1.1, 4.7, 0.1, 9.3, 0.4};
	Vector<Real>	vecSol(rhs.size());

	LUDecompositionSolver<Real> luSolver(mat);

	luSolver.Solve(rhs, vecSol);

	Vector<Real>    res_rhs = mat * vecSol;
~~~

***Functions - from function pointer***
Introduction - par recenica, sve prezentirati odjednom
~~~ c++
    double Demo_Function_TestFunc(double x) 
    { 
        return sin(x)*(1.0 + 0.5*x*x); 
    }
    ...
    // creating a function object from an already existing (standalone) function
    RealFunction f1(Demo_Function_TestFunc);

    // or creating a function object directly with lambda
    RealFunction f2{[](double x) { return sin(x)*(1.0 + 0.5*x*x); } };
~~~

Functions - from class member functions, using std::function wrapper
Introduction - par recenica
~~~ c++
    // TODO - functions example
~~~

*** Working with functions - derivation, integration***
Introduction - par recenica
~~~ c++
    // TODO - functions example
~~~

***Parametric curves - differential geometry***
~~~ c++
    // creating curve directly with lambda
    // TODO - dodati 2D i 2D polar
    ParametricCurve<3>        test_curve1( [](double t) -> VectorN<Real, 3> { return VectorN<Real, 3>{t, t*t, t*t*t}; } );
    
    // using predefined curve
    Curves::HelixCurve        helix(2.0, 2.0);
    
    // using curve from TestData
    const ParametricCurve<3> &test_curve = TestBeds::ParametricCurvesTestBed::_listCurves[0]._curve;

    double t = 0.5;
    auto tangent   = DiffGeometry::getTangent(test_curve, t);
    auto unit_tang = DiffGeometry::getTangentUnit(test_curve, t);
    auto normal    = DiffGeometry::getNormal(test_curve, t);
    auto unit_norm = DiffGeometry::getNormalUnit(test_curve, t);
    auto binormal  = VectorProd(Vector3Cartesian(unit_tang), Vector3Cartesian(unit_norm));
    
    auto curv_vec   = DiffGeometry::getCurvatureVector(test_curve, t);
    auto curvature  = DiffGeometry::getCurvature(test_curve, t);
    auto curvature3 = DiffGeometry::getCurvature3(test_curve, t);
~~~
vizualizirati neku krivulju, i u jednoj točki vizualizirati (World view) 3 vektora tangente, normale i binormale, te vektore zakrivljenosti

***Solving ODE system***
u primjeru definirati Van Der Polov oscilator, i riješiti ga
~~~ c++
    auto sys0 = TestBeds::ODESystemTestBed::getODESystem(1);

    const double atol=1.0e-3, rtol=atol;
    const double h1=0.01, hmin=0.0;
    const double x1=0.0, x2=2.0;
    
    Vector<Real> ystart01{2.0, 0.0};
    Output out(20);            

    ODESystemSolver<StepperDopr853> ode_solver01(sys0, atol, rtol, out);
    ODESystemSolution               sol01 = ode_solver01.integrate(ystart01, 0.0, 2.0, h1, hmin);

    std::cout << "x values:\n";  sol01.xval.Print(std::cout, 6, 3); std::cout << std::endl;
    std::cout << "y values: - "; sol01.yval.Print(std::cout, 6, 3);
~~~
prikazati vizualizaciju pojave chaosa za odabrane parametre

***Fields and field operations - grad, div, curl, Laplacian***
~~~ c++
    // TODO - fields example
~~~

***Analyzer - function, curve, ?***
~~~ c++
    // TODO - analyzer example
~~~

**Testing precision**

In applying any kind of numerical procedure on computers, observing precision is of paramount importance.
Following sections describe the precision of implemented algorithms.

[Intro](/docs/testing_precision/testing_precision_intro.md) - introduction

- derivation precision
- integration precision
- interpolation precision
- vector field operations precision

**Testing speed**
// TODO - testing speed
- of paramount essence, if you are using C++ (if not, you would be using Python, numpy, and scipy)

**LICENSING**
- Code is given as it is, without any warranty. Use it at your own risk.
- STRICTLY NON-COMMERCIAL USE ONLY!
- Unfortunately, also unavailable for Open Source project, due to restrictive Numerical Recipes license (for which code I have only personal license).
- So basically, it is for personal, educational and research use only.


