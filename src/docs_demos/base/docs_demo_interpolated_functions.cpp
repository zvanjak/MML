#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/InterpolatedFunction.h"
#endif

using namespace MML;

void Docs_Demo_Linear_Interpolation()
{
    std::cout << "=== Linear Interpolation ===" << std::endl << std::endl;
    
    Vector<Real> x({ 0.0, 1.0, 2.0, 3.0, 4.0 });
    Vector<Real> y({ 0.0, 1.0, 0.5, 2.0, 1.5 });
    
    LinearInterpRealFunc f(x, y, false);  // no extrapolation
    
    std::cout << "Data points: (0,0), (1,1), (2,0.5), (3,2), (4,1.5)" << std::endl;
    std::cout << "Interpolated values:" << std::endl;
    std::cout << "  f(0.5) = " << f(0.5) << std::endl;
    std::cout << "  f(1.5) = " << f(1.5) << std::endl;
    std::cout << "  f(2.5) = " << f(2.5) << std::endl << std::endl;
}

void Docs_Demo_Polynomial_Interpolation()
{
    std::cout << "=== Polynomial Interpolation ===" << std::endl << std::endl;
    
    Vector<Real> x({ 0.0, 1.0, 2.0, 3.0, 4.0, 5.0 });
    Vector<Real> y({ 0.0, 0.8, 0.9, 0.1, -0.8, -1.0 });
    
    PolynomInterpRealFunc f(x, y, 3);  // local cubic (order 3)
    
    std::cout << "Data points from smooth function" << std::endl;
    Real val = f(2.5);
    Real err = f.getLastErrorEst();
    std::cout << "  f(2.5) = " << val << ", error estimate: " << err << std::endl << std::endl;
}

void Docs_Demo_Spline_Interpolation()
{
    std::cout << "=== Cubic Spline Interpolation ===" << std::endl << std::endl;
    
    Vector<Real> x({ 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 });
    Vector<Real> y({ 0.0, 0.8, 0.9, 0.1, -0.8, -1.0, -0.5 });
    
    // Natural spline
    SplineInterpRealFunc f(x, y);
    
    std::cout << "Natural cubic spline:" << std::endl;
    std::cout << "  f(2.5) = " << f(2.5) << std::endl;
    std::cout << "  f'(2.5) = " << f.Derivative(2.5) << std::endl;
    std::cout << "  f''(2.5) = " << f.SecondDerivative(2.5) << std::endl;
    std::cout << "  Integral [1,4] = " << f.Integrate(1.0, 4.0) << std::endl << std::endl;
}

void Docs_Demo_Barycentric_Interpolation()
{
    std::cout << "=== Barycentric Rational Interpolation ===" << std::endl << std::endl;
    
    Vector<Real> x({ 0.0, 1.0, 2.0, 3.0, 4.0, 5.0 });
    Vector<Real> y({ 1.0, 1.5, 1.2, 2.0, 1.8, 2.5 });
    
    BarycentricRationalInterp f(x, y, 3);  // d=3 (default)
    
    std::cout << "Pole-free rational interpolation:" << std::endl;
    std::cout << "  f(2.5) = " << f(2.5) << " (guaranteed no pole)" << std::endl;
    
    // Verify exact at data points
    std::cout << "  Exact at data: f(x[2]) = " << f(x[2]) << ", y[2] = " << y[2] << std::endl << std::endl;
}

void Docs_Demo_Bilinear_2D()
{
    std::cout << "=== Bilinear 2D Interpolation ===" << std::endl << std::endl;
    
    // Create a 5x5 grid
    Vector<Real> x({ 0.0, 1.0, 2.0, 3.0, 4.0 });
    Vector<Real> y({ 0.0, 1.0, 2.0, 3.0, 4.0 });
    
    // z values: z(i,j) = x[i] + y[j]
    Matrix<Real> z(5, 5);
    for (int i = 0; i < 5; ++i)
        for (int j = 0; j < 5; ++j)
            z(i, j) = x[i] + y[j];
    
    BilinearInterp2D f(x, y, z);
    
    std::cout << "2D grid with z = x + y:" << std::endl;
    std::cout << "  f(1.5, 2.5) = " << f(1.5, 2.5) << " (expected: 4.0)" << std::endl;
    std::cout << "  f(2.0, 3.0) = " << f(2.0, 3.0) << " (exact at grid)" << std::endl << std::endl;
}

void Docs_Demo_Parametric_Curve()
{
    std::cout << "=== Parametric Curve Interpolation ===" << std::endl << std::endl;
    
    Matrix<Real> points(5, 3, {
        0.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        1.0, 1.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0
    });
    
    SplineInterpParametricCurve<3> curve(points, false);
    
    std::cout << "3D spline curve through 5 points:" << std::endl;
    VectorN<Real, 3> pos = curve(0.5);
    std::cout << "  Position at t=0.5: (" << pos[0] << ", " << pos[1] << ", " << pos[2] << ")" << std::endl;
    
    pos = curve(0.0);
    std::cout << "  Position at t=0: (" << pos[0] << ", " << pos[1] << ", " << pos[2] << ")" << std::endl;
    
    pos = curve(1.0);
    std::cout << "  Position at t=1: (" << pos[0] << ", " << pos[1] << ", " << pos[2] << ")" << std::endl << std::endl;
}

void Docs_Demo_Interpolated_functions()
{
    std::cout << "***********************************************************" << std::endl;
    std::cout << "*****          Interpolated Functions Demo            *****" << std::endl;
    std::cout << "***********************************************************" << std::endl << std::endl;
    
    Docs_Demo_Linear_Interpolation();
    Docs_Demo_Polynomial_Interpolation();
    Docs_Demo_Spline_Interpolation();
    Docs_Demo_Barycentric_Interpolation();
    Docs_Demo_Bilinear_2D();
    Docs_Demo_Parametric_Curve();
}