# MML - Minimal Math Library
All your basic (numerical) math needs, contained in a single-header file.

## MML Vision
- For a C++ developer, on Windows, Mac or Linux
- Who needs a math library to perform simple (and not so simple) numerical calculations 
- The Minimal Math Library is a general purpose, pythonesque in its focus on simplicity of use single-header C++ library of classes and functions
- That is trivial to use in any kind of project, is C++ 20 cross-platform compatible, and comes with a rich set of functionalities for working with vectors, matrices, tensors, linear systems, real, scalar and vector functions, polynoms, differential equations, coordinate systems and transformations, with algorithms for derivation, integration, root finding, interpolation, optimization, statistics, and more.

## Basic facts
- As of now, and for foreseable future, this is unfortunately strictly for personal, research and educational use only (see Licensing at the end)
- Single-header  - get MML.h, include in your project and go 
- Cross-platform - tested with Visual Studio Code on Windows (MSVC++, clang), Mac (gcc, clang), Linux (gcc)
- C++20 standard - but can easily be adapted to C++17, C++14
- Pythonesque in its focus on simplicity of use, and focused on faithful modeling of mathematical entities (while trying as much as possible to retain C++ computational efficiency)
- Currently, visualizators are available only on Windows platform

**Is there really a need for another C++ math library?**
- Main benefit, and the reason I did all this is attempt at completeness, and simplicity of use (also, tensors ❤️).
- Not hard focused on efficiency so when going gets rough - Boost, Eigen, Armadillo, etc.

## Organization

Library is organized in three main groups (you could call them layers), with addition of visualizers as external tools.

**Base**

Basic math types. These are the building blocks of the library, sitting at the lowest layer, depending only on standard library headers (and possibly Vector and Matrix class), and are used in all other parts of the library. 
- [Algebra](/docs/base/Algebra.md) - groups, permutation group (big TODO!)
- [Vectors](/docs/base/Vectors.md) - Vector, VectorN<int N>, Vector(2)(3)Cartesian, Vector2Polar, Vector3Spherical, Vector3Cylindrical
- [Matrices](/docs/base/Matrices.md) - Matrix, MatrixNM<int N, int M>, MatrixSym, MatrixTridiag, MatrixBandDiag
- [Tensors](/docs/base/Tensors.md) - Tensor2<int Dim>, Tensor3<int Dim>, Tensor4<int Dim>, Tensor5<int Dim> in N dimensions
- [Polynoms](/docs/base/Polynoms.md) - general Polynom class (tested for Real, Complex and Matrix as field type)
- [Geometry](/docs/base/Geometry.md) - pure geometry: points, triangles, bodies
- [2D & 3D geometry](/docs/base/Geometry_2D_3D.md) - analytic geometry in 2D and 3D
- [Vector spaces](/docs/base/Vector_spaces.md) - vector space, normed vector space, metric (Hilbert) space (still much to do here!)
- [Functionals, operators, quadratic forms](/docs/base/Operators.md)  - linear functional, quadratic form, linear operator (much to do here!)
- [Base utils](/docs/base/BaseUtils.md) - general utilities including matrix helper (IsOrthogonal, IsUnitary, IsHermitian)

**Core**

Core mathematical objects and operations of the library, depending on Base types, and used by higher algorithms.
Function objects, and different algorithms for working with them are the heart of this layer.
- [Linear alg. equations solvers](/docs/core/Linear_equations_solvers.md) - GJ, LU, QR, SVD, Cholesky
- [Functions](/docs/core/Function_types.md) - IRealFunction, IScalarFunction<int N>, IVectorFunction<int N>, IParametricCurve<N>, IParametricSurface<N>, ITensorField<N>
- [Standard functions](/docs/core/Functions.md) - definition of available standard functions
- [Interpolated functions](/docs/core/Interpolated_functions.md) - linear, polynomial, rational polynomial, spline interpolations
- [Dirac delta function](/docs/core/Dirac_delta_function.md)- predefined distributions for representing Dirac delta function
- [Numerical derivation](/docs/core/Derivation.md) - orders 1, 2, 4, 6, 8 for IRealFunction, IScalarFunction, IVectorFunction, IParametricCurve, IParametricSurface, ITensorField
- [Numerical integration](/docs/core/Integration.md) - Trapezoidal, Simpson, Romberg basic integration algorithms
- [Multidim integration](/docs/core/Multidim_integration.md) - calculating 2D and 3D (cartesian) integrals
- [Curves & Surfaces](/docs/core/Curves_and_surfaces.md) - predefined curves and surfaces
- [Fields](/docs/core/Fields.md) - predefined example fields
- [Field operations](/docs/core/Vector_field_operations.md) - grad, div, curl, Laplacian in general, Cartesian, cylindrical and spherical coordinates
- [Metric tensor](/docs/core/Metric_tensor.md) - predefined metric tensors in General, Cartesian, Cylindrical and Spherical coordinates
- [Coordinate transformations](/docs/core/Coordinate_transformations.md) - General, Cartesian, Cylindrical, Spherical  
- [ODE system](/docs/core/ODE_system.md) - represents a dynamical system of ordinary differential equations
- [Function spaces](/docs/core/Function_spaces.md) - Hermitian, Legendre, Laguerre, Chebyshev, Fourier spaces
- [Core utils](/docs/core/CoreUtils.md) - general core utilities 

**Algorithms**

Algorithms for solving mathematical problems. These are the algorithms of the library, depending on Base and Core types.
- [Eigen solvers](/docs/algorithms/Eigen_solvers.md) - solving eigenvalue problems for symmetric and non-symmetric real matrices
- [Path integration](/docs/algorithms/Path_integration.md) - calculating line and work integrals
- [ODE system solvers](/docs/algorithms/Differential_equations_solvers.md) - solvers for systems of ordinary differential equations
- [Differential geometry](/docs/algorithms/Differential_geometry.md) - for curves only, so far
- [Root finding](/docs/algorithms/Root_finding.md) - different root finding algorithms (bracketing, Newton-Raphson)
- [Statistics](/docs/algorithms/Statistics.md) - basics - avg, std, var, cov, corr
- [Function analyzer](/docs/algorithms/Function_analyzer.md) - analyzing functions, finding roots, extrema, inflection points, etc.
- [Fourier transformation](/docs/algorithms/Fourier_transformation.md) - TODO

**Visualizers**

Set of external tools, used for visualization of different types of functions. They are using as input files created by Serialize members of respective function types.
- [Real function visualizer](/docs/visualizers/RealFunction_visualizer.md) - visualizes real function in 2D
- [Parametric curve visualizer](/docs/visualizers/ParametricCurve_visualizer.md) - visualizes parametric curve in 3D
- [Surface visualizer](/docs/visualizers/Surface_visualizer.md) - visualizes cartesian 2D surface in 3D
- [Vector field visualizer](/docs/visualizers/VectorField_visualizer.md) - visualizes vector field in 3D

Example visualizations:

![My Image](docs/images/readme_visualizators_example.png)

## Simple example of real use

I decided to try something with electromagnetism, because it nicely intersects special relativity (which I need to figure out first, before going to general theory, main objective of my research) with theory of tensors and differential forms, also my great interest, and I have chosen to calculate and visualize EM field of three infinite line currents.
~~~ c++
// simple structure modeling line current
struct LineCurrent {
    double _currentI;
    Line3D _line;
};
// class representing calculated EM vector field
class InfiniteLineCurrent_Field_B : public IVectorFunction<3>
{
    std::vector<LineCurrent> _lines;
public:
    void AddLine(double currentI, const Line3D &line) { _lines.push_back({currentI, line}); }

    VectorN<Real, 3> operator()(const VectorN<Real, 3> &x) const override  {
        VectorN<Real, 3>  ret;
        Point3Cartesian   pos(x[0], x[1], x[2]);

        // calculating contribution to B field of each line at given point
        for(int i=0; i<_lines.size(); i++) {
            Point3Cartesian  nearest_point = _lines[i]._line.NearestPoint(pos);
            Vector3Cartesian vec_to_pos (nearest_point, pos);
            Vector3Cartesian fieldDirection = VectorProd(_lines[i]._line.Direction(), vec_to_pos);

            double B_magnitude = _lines[i]._currentI / (2 * Constants::PI * pos.Dist(nearest_point));

            ret = ret + B_magnitude * fieldDirection.GetAsUnitVector();
        }
        return ret;
    }
};

// main program, with setup, serialization and visualization
InfiniteLineCurrent_Field_B   EM_field;
EM_field.AddLine(300.0, Line3D(Point3Cartesian(120, 50, -50), Vector3Cartesian(0, 1, 1)));
EM_field.AddLine(200.0, Line3D(Point3Cartesian(-150, 100, 0), Vector3Cartesian(0, 0, 1)));
EM_field.AddLine(200.0, Line3D(Point3Cartesian(20, -100, 00), Vector3Cartesian(1, 0, 0)));

EM_field.Serialize3DCartesian(-300.0, 300.0, 30, -300.0, 300.0, 30, -300.0, 300.0, 30, "EM_field.txt", 3.0);
std::system("..\\..\\tools\\visualizers\\vector_field_visualizer\\MML_VectorFieldVisualizer.exe EM_field.txt");
~~~
![My Image](docs/images/readme_MainExample_EM_field.png)

## More real use examples
Before basic introductory examples, couple of real examples what it can be used for. With important note that, alas, all of them are still on ToDo list, so it is actually a plan 🙄
- at elevation 45deg, ball is fired with speed 10, 100, 1000, 10e5, 10e7 m/s, where it will be in 1 hour? [link](/docs/examples/Example1_kosi_hitac.md)
- collision calculator, 2D and 3D - [link](/docs/examples/Example2_collision_calculator.md)
- calculating tensor of inertia - [link](/docs/examples/Example3_tensor_of_inertia.md)
- investigating gravity field - [link](/docs/examples/Example4_gravity_field_investigations.md)
- Voyager travels through Solar system - [link](/docs/examples/Example5_Voyager_travels.md)
- electric charge distribution in solid body - [link](/docs/examples/Example6_electric_charge_distribution.md)

## Introductory examples

In the following sections, some basic examples of using the library are given.

**Vectors, matrices**

Examples of basic vector and matrix operations
~~~ c++
Vector<double>  vec1({ 1.5, -2.0, 0.5 }), vec2({ 1.0, 1.0, -3.0 }); 
VectorComplex   vec_cmplx1({ Complex(1,1), Complex(-1,2) });
VectorComplex   vec_cmplx2({ Complex(1,1), Complex(-1,2), Complex(2.5, -1.5) });

Matrix<double>  mat_3x3{ 3, 3, { 1.0,  2.0, -1.0, 
                                -1.0,  5.0,  6.0, 
                                 3.0, -2.0,  1.0 }};  
MatrixComplex   mat_cmplx(2,2, { Complex(0.5,1), Complex(-1,2), 
                                 Complex(-1,-2), Complex(-2,2) });
MatrixComplex   mat_cmplx2(2,3, { Complex(1,2),    Complex(-1,1), Complex(1.5,-2), 
                                  Complex(2,-0.5), Complex(3,-2), Complex(-1,1) });
Matrix<double>  unit_mat3 = MML::Matrix<Real>::GetUnitMatrix(3);

Vector<double> v_real  = 2.0 * (vec1 + vec2) * mat_3x3 / vec1.NormL2();
VectorComplex  v_cmplx = vec_cmplx1 * mat_cmplx / Complex(1.5, -1.5) / 2.0;

std::cout << "v_real  = " << v_real << std::endl;
std::cout << "v_cmplx = " << v_cmplx << std::endl;

// combining real and complex vectors and matrices requires special functions
VectorComplex cvec2 = MatrixUtils::MulVecMat(vec_cmplx2, mat_3x3);
MatrixComplex cmat2 = MatrixUtils::MulMat(mat_cmplx2, mat_3x3);

std::cout << "Matrix mat_3x3   = " << mat_3x3 << std::endl;
std::cout << "Matrix mat_cmplx = " << mat_cmplx << std::endl;

std::cout << "IsOrthogonal(mat_3x3)  = " << MatrixUtils::IsOrthogonal(mat_3x3) << std::endl;
std::cout << "IsHermitian(mat_cmplx) = " << MatrixUtils::IsHermitian(mat_cmplx) << std::endl;
std::cout << "IsUnitary(mat_cmplx)   = " << MatrixUtils::IsUnitary(mat_cmplx) << std::endl;

/* OUTPUT
    v_real  = [   -3.137858162,     3.922322703,    -8.629109946]
    v_cmplx = [        (0.5,1), (-0,-1.666666667)]
    Matrix mat_3x3   = Rows: 3 Cols: 3
    [          1,          2,         -1,  ]
    [         -1,          5,          6,  ]
    [          3,         -2,          1,  ]

    Matrix mat_cmplx = Rows: 2 Cols: 2
    [    (0.5,1),     (-1,2),  ]
    [    (-1,-2),     (-2,2),  ]

    IsOrthogonal(mat_3x3)  = 0
    IsHermitian(mat_cmplx) = 1
    IsUnitary(mat_cmplx)   = 0
*/
~~~

**Solving linear systems of equations and calculating eigenvalues**

How to solve linear systems of equations, and calculate eigenvalues
~~~ c++
Matrix<Real>    mat{5, 5, { 0.2,  4.1, -2.1, -7.4,  1.6,
                            1.6,  1.5, -1.1,  0.7,  5.0,
                           -3.8, -8.0,  9.6, -5.4, -7.8,
                            4.6, -8.2,  8.4,  0.4,  8.0,
                           -2.6,  2.9,  0.1, -9.6, -2.7 } };
Vector<Real> 	rhs{1.1, 4.7, 0.1, 9.3, 0.4};

LUDecompositionSolver<Real> luSolver(mat);

Vector<Real>	vecSol = luSolver.Solve(rhs);

std::cout << "Solution   = " << vecSol << std::endl;
std::cout << "Right side = "; rhs.Print(std::cout,8,4); std::cout << std::endl;
std::cout << "Mat * sol  = "; (mat * vecSol).Print(std::cout,8,4); std::cout << std::endl;

Matrix<Real>  matcopy(mat);

EigenSolver   eigenSolver(matcopy, true, false);

std::cout << "\nNum real eigenvalues    : " << eigenSolver.getNumReal();
std::cout << "\nNum complex eigenvalues : " << eigenSolver.getNumComplex() << "\n";

std::cout << "\nEigenvalues : "; eigenSolver.getEigenvalues().Print(std::cout,10,5); 
std::cout << "\nReal        : "; eigenSolver.getRealEigenvalues().Print(std::cout,15,10); 
std::cout << "\nComplex     : "; eigenSolver.getComplexEigenvalues().Print(std::cout,15,10); 

/* OUTPUT
    Solution   = [   -5.568500786,    -5.944693206,    -5.007620645,    -1.393638021,     3.598760994]
    Right side = [     1.1,      4.7,      0.1,      9.3,      0.4]
    Mat * sol  = [     1.1,      4.7,      0.1,      9.3,      0.4]

    Num real eigenvalues    : 3
    Num complex eigenvalues : 2

    Eigenvalues : [(12.974,0), (0.99944,0), (-0.033184,0), (-2.4701,12.994), (-2.4701,-12.994)]
    Real        : [    12.97392154,    0.9994371124,  -0.03318390189]
    Complex     : [(-2.470087376,12.99433106), (-2.470087376,-12.99433106)]
*/
~~~

**Defining functions**

Four (main) possibilities for defining/creating functions of various types
~~~ c++
// CASE 1 - standalone function providing calculation of a function
double Readme_functions_TestFunc(double x) { 
    return sin(x)*(1.0 + 0.5*x*x); 
}

// creating a function object from an already existing (standalone) function
RealFunction f1(Readme_functions_TestFunc);

// CASE 2 - create it directly with lambda
RealFunction f2{[](double x) { return sin(x)*(1.0 + 0.5*x*x); } };

// creating directly different types of functions
ScalarFunction<3>       funcScalar([](const VectorN<Real, 3> &x) { return x[0]; });
VectorFunction<3>       funcVector([](const VectorN<Real, 3> &x) { return VectorN<Real, 3>{0, x[0] * x[1], 0}; });
VectorFunctionNM<2, 3>  funcVectorNM([](const VectorN<Real, 2> &x) { return VectorN<Real, 3>{0, x[0] * x[1], 0}; });
ParametricCurve<3>      paramCurve([](double x) { return VectorN<Real, 3>{x, 2 * x, 3 * x}; });
ParametricSurface<3>    paramSurface([](double x, double y) { return VectorN<Real, 3>{x * y, 2 * x * y, 3 * x}; });  

// CASE 3 - class you CAN change has member function that does the calculation
// Option 1 - define operator() for your class and create RealFunctionFromStdFunc
class ClassProvidingFuncToDerive {
    public:
        double operator()(double x ) const { 
            return 1.0;     /* calculation using member variables */ 
        }
};

ClassProvidingFuncToDerive   obj1;
RealFunctionFromStdFunc f1(std::function<double(double)>{obj1});

// Option 2 - make your class inherit IRealFunction interface, and use the object itself as RealFunction
class ClassProvidingFuncToDerive2 : public IRealFunction {
    public:
        double operator()(double x ) const { 
            return 1.0;     /* calculation using member variables */ 
        }
};

ClassProvidingFuncToDerive2   f2;       // usable RealFunction object (can be derived, integrated, ...)

// CASE 4 - class you CAN'T change has member function that does the calculation
class BigComplexClassYouCantChange {
    // has data for calculating function you want to do somethign with
};

// Create a helper wrapper class, inherit it from IRealFunction and use it as RealFunction
class BigComplexFunc2 : public IRealFunction {
    const BigComplexClassYouCantChange &_ref;
public:
    BigComplexFunc2(const BigComplexClassYouCantChange &bigClass) : _ref(bigClass) { }

    double operator()(double x ) const {
        return 1.0;     /* calculation using _ref */ 
    }
};

BigComplexClassYouCantChange bigObj;
BigComplexFunc2    f1(bigObj);          // usable RealFunction object (can be derived, integrated, ...)
~~~

**Creating interpolated functions**

If there is no explicit function, but we have data, there are various ways to create interpolated function
~~~ c++
const int NumInterpPnt = 5;
const Real x1 = 0, x2 = 10.0;

// we will use this as test func
RealFunction test_func{[](double x) { return sin(x)*(1.0 + 0.5*x*x); } };

// and using our helper, available for all real functions, create data for interpolation
Vector<double> x_val(NumInterpPnt), y_val(NumInterpPnt);
test_func.GetValues(x1, x2, NumInterpPnt, x_val, y_val);

// these are the ways we can interpolate Real function
LinearInterpRealFunc    f_linear(x_val, y_val);
PolynomInterpRealFunc   f_polynom(x_val, y_val, 3);
SplineInterpRealFunc    f_spline(x_val, y_val);
BaryRatInterpRealFunc   f_baryrat(x_val, y_val, 3);

test_func.SerializeEquallySpacedDetailed(x1, x2, 100, "..\\..\\results\\readme_interp_test_func.txt");
f_linear.SerializeEquallySpacedDetailed(x1, x2, 500, "..\\..\\results\\readme_interp_linear_5_pnt.txt");
f_polynom.SerializeEquallySpacedDetailed(x1, x2, 100, "..\\..\\results\\readme_interp_polynom_5_pnt.txt");
f_spline.SerializeEquallySpacedDetailed(x1, x2, 100, "..\\..\\results\\readme_interp_spline_5_pnt.txt");
f_baryrat.SerializeEquallySpacedDetailed(x1, x2, 100, "..\\..\\results\\readme_interp_baryrat_5_pnt.txt");

const char *cmd = "..\\..\\tools\\visualizers\\real_function_visualizer\\MML_RealFunctionVisualizer.exe"
                    " ..\\..\\results\\readme_interp_test_func.txt"
                    " ..\\..\\results\\readme_interp_linear_5_pnt.txt"
                    " ..\\..\\results\\readme_interp_polynom_5_pnt.txt"
                    " ..\\..\\results\\readme_interp_spline_5_pnt.txt"
                    " ..\\..\\results\\readme_interp_baryrat_5_pnt.txt";
std::system(cmd);
~~~
Visualization of test func, and interpolations for NumInterpPnt = 5, 8 and 12

![My Image](docs/images/readme_interpolated_functions.png)

**Working with functions - derivation**

Examples of various ways of calculating derivation of functions
~~~ c++
RealFunction       f1{[](double x) { return sin(x)*(1.0 + 0.5*x*x); } };

// numerical derivation of real function (available orders - 1, 2, 4, 6, 8)
double der_f1 = Derivation::NDer1(f1, 0.5);
double der_f4 = Derivation::NDer2(f1, 0.5, 1e-6);   // setting explicit step size
double err;
double der_f6 = Derivation::NDer6(f1, 0.5, &err);   // if we need error estimate    
// we can use default Derive routine (set to NDer4), but it requires error estimate
double num_der4 = Derivation::Derive(f1, 0.5, nullptr);

// second and third derivatives
double sec_der_f1   = Derivation::NSecDer2(f1, 0.5);
double third_der_f1 = Derivation::NThirdDer2(f1, 0.5);

// creating new function that is derivation of existing RealFunction
RealFuncDerived4    f1_der4(f1);        // 4th order derivation

// scalar and vector functions
ScalarFunction<3>   f2Scal([](const VectorN<Real, 3> &x) { return 1.0 / pow(x.NormL2(), 2); });
VectorFunction<3>   f3Vec([](const VectorN<Real, 3> &x) { return VectorN<Real, 3>{0, x[0] * x[1], 0}; });
VectorN<Real, 3>    der_point{1.0, 1.0, 1.0};

double der_f2               = Derivation::NDer1Partial(f2Scal, 1, der_point);
VectorN<Real, 3> der_f2_all = Derivation::NDer1PartialByAll(f2Scal, der_point);

double der_f3 = Derivation::NDer1Partial(f3Vec, 1, 1, der_point);
VectorN<Real, 3>     der_f3_by1    = Derivation::NDer2PartialByAll(f3Vec, 1, der_point);
MatrixNM<Real, 3, 3> der_f3_by_all = Derivation::NDer4PartialAllByAll(f3Vec, der_point);
~~~

**Working with functions - integration**

Examples of integration real functions, but also 2D and 3D scalar functions
~~~ c++
// numerical integration of real function
RealFunction f1{[](double x) { return sin(x)*(1.0 + 0.5*x*x); } };

double a = 0.0;
double b = 1.0;
double int_trap = IntegrateTrap(f1,a,b);
double int_simp = IntegrateSimpson(f1,a,b);
double int_romb = IntegrateRomberg(f1,a,b);
// we can use default Integrate routine (set to IntegrateSimpson), requires precision
double int_def = Integrate(f1, a, b, 1e-04);

// 2D integration of constant scalar 2D function (ie. we'll get the area of the surface)
ScalarFunction<2> f2([](const VectorN<Real, 2> &x) { return 1.0; });

// we integrate over circle with radius 2
Real val = IntegrateSurface(f2, IntegrationMethod::GAUSS10, 
                                -2, 2,              // x range
                                [](Real x) { return -sqrt(4 - x*x);},   // y range lower limit
                                [](Real x) { return sqrt(4 - x*x);});   // y range upper limit

std::cout << "Calc. area = " << val << ", exact value: 4 * PI = " << 4 * Constants::PI << std::endl;

// 3D integration of constant scalar 3D function (ie. we'll get the volume of the solid)
ScalarFunction<3> f3([](const VectorN<Real, 3> &x) { return 1.0; });

Real vol = IntegrateVolume( f3, 
                            -1, 1,              
                            [](Real x) { return -sqrt(1 - x*x);}, 
                            [](Real x) { return sqrt(1 - x*x);}, 
                            [](Real x, Real y) { return -sqrt(1 - x*x - y*y);}, 
                            [](Real x, Real y) { return sqrt(1 - x*x - y*y);});

std::cout << "Calc. vol. = " << vol << ", exact value: 4/3 * PI = " << 4.0/3.0 * Constants::PI << std::endl;

/* OUTPUT 
    Calc. area = 12.57211164, exact value: 4 * PI = 12.56637061
    Calc. vol. = 4.190703882, exact value: 4/3 * PI = 4.188790205
*/
~~~

**Solving ODE system**

Solving system of ordinary differential equations
~~~ c++
// in-place definition of Lorenz system (have to fix parameters!)
ODESystem LorenzSystem(3, [](double t, const Vector<Real>& y, Vector<Real>& dydt)
{
    double sigma = 10.0, rho = 28.0, beta = 8.0 / 3.0;
    dydt[0] = sigma * (y[1] - y[0]);

    dydt[1] = y[0] * (rho - y[2]) - y[1];
    dydt[2] = y[0] * y[1] - beta * y[2];
});    

// get it from predefined test-bed 
TestBeds::LorenzSystemODE alsoLorenzSystem(10.0, 28.0, 8.0/3.0);

const double atol=1.0e-3, rtol=atol, h1=0.01, hmin=0.0;
double x1=0.0, x2=50.0;

Vector<Real> init_cond({2.0, 1.0, 1.0});
Output out0(10000);

ODESystemSolver<StepperDopr5> ode_solver0(LorenzSystem,atol,rtol, out0);
ODESystemSolution             sol0 = ode_solver0.integrate(init_cond, x1, x2, h1, hmin);

sol0.Serialize("demo5_lorenz_system.txt", "Lorenz system");
auto ret2 = std::system("..\\..\\tools\\visualizers\\real_function_visualizer\\MML_RealFunctionVisualizer.exe demo5_lorenz_system.txt");

auto curve = sol0.getSolutionAsParametricCurve<3>();
sol0.SerializeAsParametricCurve3D("demo5_lorenz_system_as_parametric_curve.txt", "Lorenz system as parametric curve");
std::system("..\\..\\tools\\visualizers\\parametric_curve_visualizer\\MML_ParametricCurveVisualizer.exe demo5_lorenz_system_as_parametric_curve.txt");
~~~
Resulting images

![My Image](docs/images/readme_ODEsolvers.png)

**Fields and field operations - grad, div, curl, Laplacian**

Using as example inverse radial field, with its potential and force field, demonstrate calculation of field gradient, divergence, curl and Laplacian
Calculations are performed in Cartesian and spherical coordinates, along circle in XZ-plane, and covariant vector transformation is also demonstrated
~~~ c++
// Setting up fields and creating scalar potential and vector force field from predefined functions
static ScalarFunction<3> pot_cart_exact([](const VectorN<Real, 3> &x_cart) -> Real   { return -InverseRadialPotentialFieldCart(x_cart); });
static ScalarFunction<3> pot_spher_exact([](const VectorN<Real, 3> &x_spher) -> Real { return -InverseRadialPotentialFieldSpher(x_spher); });
static VectorFunction<3> force_field_cart_exact([](const VectorN<Real, 3> &x_cart)   { return InverseRadialPotentialForceFieldCart(x_cart); });
static VectorFunction<3> force_field_spher_exact([](const VectorN<Real, 3> &x_spher) { return InverseRadialPotentialForceFieldSph(x_spher); });

// if we have only potential, we can numerical calculate force field from it
static VectorFunction<3> force_field_cart_from_grad{ [](const VectorN<Real, 3> &x_cart)  { return -ScalarFieldOperations::GradientCart<3>(pot_cart_exact, x_cart); } };  
static VectorFunction<3> force_field_spher_from_grad([](const VectorN<Real, 3> &x_spher) { return -ScalarFieldOperations::GradientSpher(pot_spher_exact,x_spher); });

// if we have potential in one coord. system, and we need force field in another, we can calculate it 
// by first calculating force field in the same coord. system as potential, and then transforming it covariantly to desired coordinates
VectorFunction<3> force_field_cart_from_spher_pot{ [](const VectorN<Real, 3> &x_cart)
{
    auto force_spher     = -ScalarFieldOperations::GradientSpher(pot_spher_exact, CoordTransfCartToSpher.transf(x_cart));
    VectorN<Real, 3> ret = CoordTransfSpherToCart.transfVecCovariant(force_spher, x_cart);
    return ret;        // returned force field vector is now in Cartesian coordinates
} };
VectorFunction<3> force_field_spher_from_cart_pot{ [](const VectorN<Real, 3> &x_spher)
{
    auto force_cart      = -ScalarFieldOperations::GradientCart<3>(pot_cart_exact, CoordTransfSpherToCart.transf(x_spher));
    VectorN<Real, 3> ret = CoordTransfCartToSpher.transfVecCovariant(force_cart, x_spher);
    return ret;        // returned force field vector is now in Spherical coordinates
} };
~~~
~~~ c++
// calculating potential and force around circle
Curves::Circle3DXZ circle(10.0);

// calculating field Gradient in Cart. and Spherical coordinates around circle
std::cout << "         Position          Pot.(Cart)    Pot.(Spher)        Force exact (Cart. vector)              Force exact (Sph. vector)              Force num.grad (Cart.vector)        Force num.grad (Sph.vector)" << std::endl;
for(double t=0.0; t<2*Constants::PI; t+=1)
{
    Vector3Cartesian   pos_cart = circle(t);
    Vector3Spherical   pos_spher = CoordTransfCartToSpher.transf(pos_cart);

    double pot_cart  = pot_cart_exact(pos_cart);
    double pot_spher = pot_spher_exact(pos_spher);

    Vector3Cartesian   force_cart_exact    = force_field_cart_exact(pos_cart);
    Vector3Spherical   force_spher_exact   = force_field_spher_exact(pos_spher);
    Vector3Cartesian   force_cart_numgrad  = force_field_cart_from_grad(pos_cart);
    Vector3Spherical   force_spher_numgrad = force_field_spher_from_grad(pos_spher);

    std::cout << pos_cart.to_string(6,3) << "       " << pot_cart << "           "  << pot_spher << "    "
            << "  "  << force_cart_exact.to_string(10,4) 
            << "  "  << force_spher_exact.to_string(10,4) 
            << "  "  << force_cart_numgrad.to_string(10,4) 
            << "  "  << force_spher_numgrad.to_string(10,4) 
            << std::endl;
}
~~~
~~~ c++
// calculating field Divergence of four defined vector fields around circle
std::cout << "            Position               Cart.Div of               Cart.Div. of         Sph.Div of              Sph.Div of   " << std::endl;
std::cout << "                                  exact Cart.field          num.grad.field     exact Spher. field    num.grad.sph.field" << std::endl;
for(double t=0.0; t<2*Constants::PI; t+=1.0)
{
    Vector3Cartesian   pos_cart  = circle(t);
    Vector3Spherical   pos_spher = CoordTransfCartToSpher.transf(pos_cart);

    auto div_cart_exact  = VectorFieldOperations::DivCart<3>(force_field_cart_exact, pos_cart);
    auto div_spher_exact = VectorFieldOperations::DivSpher(force_field_spher_exact, pos_spher);
    
    auto div_grad_cart_numgrad  = VectorFieldOperations::DivCart<3>(force_field_cart_from_grad, pos_cart);
    auto div_grad_spher_numgrad = VectorFieldOperations::DivSpher(force_field_spher_from_grad, pos_spher);
    
    std::cout << pos_cart.to_string(7,3) << " " << std::setw(20) << div_cart_exact << "      " << std::setw(20) << div_grad_cart_numgrad << "      " 
                                            << div_spher_exact << "      " << div_grad_spher_numgrad << std::endl;
}
~~~
~~~ c++
// calculating vector field Curl and field Laplacian
std::cout << "            Position                    Cartesian curl                        Spherical curl                   Cartesian Laplacian     Spherical Laplacian" << std::endl;
for(double t=0.0; t<2*Constants::PI; t+=1)
{
    Vector3Cartesian   pos = circle(t);
    Vector3Spherical   pos_spher = CoordTransfCartToSpher.transf(pos);

    Vector3Cartesian   curl_cart  = VectorFieldOperations::CurlCart(force_field_cart_exact, pos);
    Vector3Spherical   curl_spher = VectorFieldOperations::CurlSpher(force_field_spher_exact, pos_spher);

    double   lapl_cart  = ScalarFieldOperations::LaplacianCart(pot_cart_exact, pos);
    double   lapl_spher = ScalarFieldOperations::LaplacianSpher(pot_spher_exact, pos_spher);

    std::cout << pos.to_string(7,3) << " = " << curl_cart.to_string(10,4) 
                                    << "  "  << curl_spher.to_string(10,4)
                                    << "        " << std::setw(16) << lapl_cart
                                    << "        " << std::setw(16) << lapl_spher
                                    << std::endl;
}
/* OUTPUT
         Position          Pot.(Cart)    Pot.(Spher)        Force exact (Cart. vector)              Force exact (Sph. vector)              Force num.grad (Cart.vector)        Force num.grad (Sph.vector)
[    10,      0,      0]       -0.1           -0.1      [     -0.01,         -0,         -0]  [     -0.01,          0,          0]  [     -0.01,         -0,         -0]  [     -0.01,         -0,         -0]
[   5.4,      0,   8.41]       -0.1           -0.1      [ -0.005403,         -0,  -0.008415]  [     -0.01,          0,          0]  [ -0.005403,         -0,  -0.008415]  [     -0.01,         -0,         -0]
[ -4.16,      0,   9.09]       -0.1           -0.1      [  0.004161,         -0,  -0.009093]  [     -0.01,          0,          0]  [  0.004161,         -0,  -0.009093]  [     -0.01,         -0,         -0]
[  -9.9,      0,   1.41]       -0.1           -0.1      [    0.0099,         -0,  -0.001411]  [     -0.01,          0,          0]  [    0.0099,         -0,  -0.001411]  [     -0.01,         -0,         -0]
[ -6.54,      0,  -7.57]       -0.1           -0.1      [  0.006536,         -0,   0.007568]  [     -0.01,          0,          0]  [  0.006536,         -0,   0.007568]  [     -0.01,         -0,         -0]
[  2.84,      0,  -9.59]       -0.1           -0.1      [ -0.002837,         -0,   0.009589]  [     -0.01,          0,          0]  [ -0.002837,         -0,   0.009589]  [     -0.01,         -0,         -0]
[   9.6,      0,  -2.79]       -0.1           -0.1      [ -0.009602,         -0,   0.002794]  [     -0.01,          0,          0]  [ -0.009602,         -0,   0.002794]  [     -0.01,         -0,         -0]
            Position               Cart.Div of               Cart.Div. of         Sph.Div of              Sph.Div of
                                  exact Cart.field          num.grad.field     exact Spher. field    num.grad.sph.field
[     10,       0,       0]     -6.505213035e-18          -9.889788207e-12      -9.728329253e-16      -2.036693314e-12
[    5.4,       0,    8.41]      1.140363845e-15          -1.403275347e-11      -9.728329253e-16      -2.036693314e-12
[  -4.16,       0,    9.09]     -2.109857428e-15          -2.632835248e-11      -9.728329253e-16      -2.036693314e-12
[   -9.9,       0,    1.41]       1.86168355e-15           7.353272722e-13      -9.728329253e-16      -2.036693314e-12
[  -6.54,       0,   -7.57]     -7.254396736e-16          -1.703953005e-11      -9.728329253e-16      -2.036693314e-12
[   2.84,       0,   -9.59]     -3.102769777e-15          -1.202754346e-11      -9.728329253e-16      -2.036693314e-12
[    9.6,       0,   -2.79]      7.181755191e-16           5.078564305e-12      -9.728329253e-16      -2.036693314e-12
            Position                    Cartesian curl                        Spherical curl                   Cartesian Laplacian     Spherical Laplacian
[     10,       0,       0] = [         0,          0,          0]  [         0,          0,          0]        -2.865156029e-14         7.464116131e-13
[    5.4,       0,    8.41] = [         0,  2.408e-15,          0]  [         0,          0,          0]        -3.885927604e-12         7.464116131e-13
[  -4.16,       0,    9.09] = [         0,  1.806e-16,          0]  [         0,          0,          0]         5.927392041e-12         7.464116131e-13
[   -9.9,       0,    1.41] = [         0, -1.715e-15,          0]  [         0,          0,          0]         5.756142091e-13         7.464116131e-13
[  -6.54,       0,   -7.57] = [         0,  6.622e-16,          0]  [         0,          0,          0]        -3.369075635e-12         7.464116131e-13
[   2.84,       0,   -9.59] = [         0,  1.625e-15,          0]  [         0,          0,          0]        -2.460620032e-12         7.464116131e-13
[    9.6,       0,   -2.79] = [         0, -2.769e-15,          0]  [         0,          0,          0]         6.892316579e-13         7.464116131e-13
*/
~~~
**Parametric curves - basic differential geometry**

Calculating tangent, normal, binormal, curvature, curvature3, curvature vector for Parametric curve.
~~~ c++
// creating curve directly with lambda
ParametricCurve<3>        test_curve1( [](double t) -> VectorN<Real, 3> { return VectorN<Real, 3>{t, t*t, t*t*t}; } );

// using predefined curve
Curves::LemniscateCurve   lemniscate;
Curves::ToroidalSpiralCurve torus(3, 2.0);

// using curve from TestData
const ParametricCurve<3> &test_curve = TestBeds::ParametricCurvesTestBed::getTestCurve("Helix")._curve;

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

## Visualizators examples

Examples how to use four available visualization tools.

**Real function visualizer**
~~~ c++
RealFunction f1{[](double x) { return sin(x) * (x-3)*(x+5) / sqrt(std::abs(2 - x)); } };

f1.SerializeEquallySpacedDetailed(-10.0, 10.0, 500, "..\\..\\results\\readme_func1.txt");
std::system("..\\..\\tools\\visualizers\\real_function_visualizer\\MML_RealFunctionVisualizer.exe ..\\..\\results\\readme_func1.txt");

RealFuncDerived4 f2(f1);        // derivation of f1 function

f2.SerializeEquallySpacedDetailed(-10.0, 10.0, 100, "..\\..\\results\\readme_func2.txt");
std::system("..\\..\\tools\\visualizers\\real_function_visualizer\\MML_RealFunctionVisualizer.exe ..\\..\\results\\readme_func2.txt");

// shown together
std::system("..\\..\\tools\\visualizers\\real_function_visualizer\\MML_RealFunctionVisualizer.exe ..\\..\\results\\readme_func1.txt ..\\..\\results\\readme_func2.txt");
~~~
Visualization:
![My Image](docs/images/readme_visualizator_real_function.png)

**Surface visualizer**

~~~ c++
// Monkey saddle surface
ScalarFunction<2> testFunc1{[](const VectorN<Real, 2> &x) { return x[0] * (x[0]*x[0] - 3 * x[1]*x[1]); } };

testFunc1.Serialize2DCartesian(-2.0, 2.0, 20, -2.0, 2.0, 20, "..\\..\\results\\readme_surface1.txt");
std::system("..\\..\\tools\\visualizers\\scalar_function_2d_visualizer\\MML_ScalarFunction2Visualizer.exe ..\\..\\results\\readme_surface1.txt");

ScalarFunction<2> testFunc2{[](const VectorN<Real, 2> &x) { return (std::abs(x[0])-10) * (std::abs(x[1])-10) * sin(x[0]) * cos(x[1]); } };

testFunc2.Serialize2DCartesian(-10.0, 10.0, 50, -10.0, 10.0, 50, "..\\..\\results\\readme_surface2.txt");
std::system("..\\..\\tools\\visualizers\\scalar_function_2d_visualizer\\MML_ScalarFunction2Visualizer.exe ..\\..\\results\\readme_surface2.txt");    
~~~
Visualization:
![My Image](docs/images/readme_visualizator_surfaces.png)

**Parametric curve visualizer**

~~~ c++
// using predefined 3D curves for visualization example
Curves::HelixCurve              helix(20.0, 2.0);
Curves::ToroidalSpiralCurve     toroid(20.0);

helix.SerializeCartesian3D(-50.0, 50.0, 1000, "readme_curve_helix.txt");
std::system("..\\..\\tools\\visualizers\\parametric_curve_visualizer\\MML_ParametricCurveVisualizer.exe readme_curve_helix.txt");

toroid.SerializeCartesian3D(0.0, 2 * Constants::PI, 5000, "readme_curve_toroid.txt");
std::system("..\\..\\tools\\visualizers\\parametric_curve_visualizer\\MML_ParametricCurveVisualizer.exe readme_curve_toroid.txt");
~~~
Visualization:
![My Image](docs/images/readme_visualizator_parametric_curve.png)

**Vector field visualizer**

~~~ c++
// defining force field of two masses as VectorFunction
VectorFunction<3> gravity_force_field{ [](const VectorN<Real, 3> &x) 
{
    const VectorN<Real, 3> x1{ 100.0, 0.0, 0.0 }, x2{ -100.0, 0.0, 0.0 };
    const Real m1 = 1000.0, m2 = 1000.0;
    const Real G = 10;
    return -G * m1 * (x - x1) / std::pow((x - x1).NormL2(), 3) - G * m2 * (x - x2) / std::pow((x - x2).NormL2(), 3);
} };

gravity_force_field.Serialize3DCartesian(-200.0, 200.0, 15, -200.0, 200.0, 15, -200.0, 200.0, 15, "readme_vector_field.txt", 5);
std::system("..\\..\\tools\\visualizers\\vector_field_visualizer\\MML_VectorFieldVisualizer.exe readme_vector_field.txt");
~~~
Visualization:
![My Image](docs/images/readme_visualizator_vector_field.png)


## Testing and precision

**Tests**

There is a set of tests for each of the above algorithms, in the /tests folder.

**Demos**

Various examples of using the library are given in the /src folder.

**Test beds**

There is a set of test beds for each of the above algorithms, in the /test_data folder.
They are used for verification of implemented algorithms, but can also be used as predefined inputs for your investigations.

- [Functions](/docs/testbeds/Functions_testbed.md) - real, scalar and vector functions test bed
- [Parametric curves & surfaces](/docs/testbeds/ParametricCurvesSurfaces_testbed.md) - parametric curves & surfaces test bed
- [Lin.alg.systems](/docs/testbeds/LinAlgSystems_testbed.md) - linear systems of equations test bed
- [ODE system](/docs/testbeds/ODESystems_testbed.md) - ODE systems test bed

**Evaluating algorithms precision and correctnes**

In applying any kind of numerical procedure on computers, observing precision is of paramount importance.
Whenever possible, implemented algorithms calculate their error estimates, and provide them as (if possible, optional) output.
Default precisions for algorithms are set in [Defaults](/docs/testing_precision/DefaultPrecisions.md) namespace, and can be changed by user.

- [Derivation precision](/docs/testing_precision/TestDerivationPrecision.md) - investigating precision of numerical derivation
- [Integration precision](/docs/testing_precision/TestIntegrationPrecision.md) - investigating precision of numerical integration
- [Interpolation precision](/docs/testing_precision/TestInterpolationPrecision.md) - investigating precision of interpolation
- [Vector fields operations precision](/docs/testing_precision/TestVectorFieldOperationsPrecision.md) - investigating precision vector field operations

## LICENSING
- Code is given as it is, without any warranty. Use it at your own risk (and convenience).
- STRICTLY NON-COMMERCIAL USE ONLY!
- Unfortunately, also unavailable for Open Source project, due to restrictive Numerical Recipes license (for which code I have only personal license).
- So basically, it is for personal, educational and research use only.
