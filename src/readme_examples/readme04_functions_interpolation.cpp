///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        readme04_functions_interpolation.cpp                                ///
///  Description: README Defining Functions & Interpolation section demo              ///
///               Demonstrates all ways to create functions and interpolation         ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///////////////////////////////////////////////////////////////////////////////////////////

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/Vector/Vector.h"
#include "base/Function.h"
#include "base/InterpolatedFunction.h"
#endif

using namespace MML;

// ==================== CASE 1: Standalone function ====================
double MyStandaloneFunction(double x) { 
    return std::sin(x) * (1.0 + 0.5*x*x); 
}

// ==================== CASE 3a: Class with operator() ====================
class FunctionFromClassOperator {
    double _amplitude;
public:
    FunctionFromClassOperator(double amp) : _amplitude(amp) {}
    double operator()(double x) const { 
        return _amplitude * std::sin(x); 
    }
};

// ==================== CASE 3b: Class inheriting IRealFunction ====================
class MyDerivedFunction : public IRealFunction {
    double _frequency;
public:
    MyDerivedFunction(double freq) : _frequency(freq) {}
    double operator()(double x) const override { 
        return std::cos(_frequency * x); 
    }
};

// ==================== CASE 4: Wrapper for external class ====================
class ExternalComplexClass {
public:
    double ComputeValue(double x) const { return std::exp(-x*x); }
};

class ExternalClassWrapper : public IRealFunction {
    const ExternalComplexClass& _ref;
public:
    ExternalClassWrapper(const ExternalComplexClass& obj) : _ref(obj) {}
    double operator()(double x) const override {
        return _ref.ComputeValue(x);
    }
};

void Readme_DefiningFunctions()
{
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                 README - Defining Functions                   ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    std::cout << "\n=== CASE 1: From standalone function ===" << std::endl;
    RealFunction f1(MyStandaloneFunction);
    std::cout << "f1(1.0) = " << f1(1.0) << std::endl;

    std::cout << "\n=== CASE 2: Direct lambda creation (most common) ===" << std::endl;
    RealFunction f2{[](double x) { return std::sin(x) * (1.0 + 0.5*x*x); }};
    std::cout << "f2(1.0) = " << f2(1.0) << std::endl;

    // Different function types with lambdas
    ScalarFunction<3> scalar([](const VectorN<Real,3>& x) { 
        return x[0]*x[0] + x[1]*x[1] + x[2]*x[2];   // r²
    });

    VectorFunction<3> vector([](const VectorN<Real,3>& x) -> VectorN<Real,3> { 
        return {x[1]*x[2], x[0]*x[2], x[0]*x[1]};   // cross-product like
    });

    ParametricCurve<3> helix([](Real t) -> VectorN<Real,3> { 
        return {std::cos(t), std::sin(t), 0.2*t}; 
    });

    std::cout << "scalar({1,2,3}) = " << scalar(VectorN<Real,3>{1,2,3}) << std::endl;
    std::cout << "vector({1,2,3}) = " << vector(VectorN<Real,3>{1,2,3}) << std::endl;
    std::cout << "helix(PI)       = " << helix(Constants::PI) << std::endl;

    std::cout << "\n=== CASE 3a: From class with operator() ===" << std::endl;
    FunctionFromClassOperator obj(2.5);  // amplitude = 2.5
    RealFunctionFromStdFunc f3(std::function<double(double)>{obj});
    std::cout << "f3(PI/2) = " << f3(Constants::PI/2) << " (expected: 2.5)" << std::endl;

    std::cout << "\n=== CASE 3b: Class inheriting IRealFunction ===" << std::endl;
    MyDerivedFunction f4(2.0);  // frequency = 2.0 → cos(2x)
    std::cout << "f4(PI/4) = " << f4(Constants::PI/4) << " (expected: 0)" << std::endl;

    std::cout << "\n=== CASE 4: Wrapper for external class ===" << std::endl;
    ExternalComplexClass externalObj;
    ExternalClassWrapper f5(externalObj);
    std::cout << "f5(0) = " << f5(0) << " (expected: 1)" << std::endl;
    std::cout << "f5(1) = " << f5(1) << " (expected: " << std::exp(-1) << ")" << std::endl;

/* Expected OUTPUT:
    === CASE 1: From standalone function ===
    f1(1.0) = 1.261465584

    === CASE 2: Direct lambda creation (most common) ===
    f2(1.0) = 1.261465584
    scalar({1,2,3}) = 14
    vector({1,2,3}) = [6, 3, 2]
    helix(PI)       = [-1, 1.22465e-16, 0.628319]

    === CASE 3a: From class with operator() ===
    f3(PI/2) = 2.5 (expected: 2.5)

    === CASE 3b: Class inheriting IRealFunction ===
    f4(PI/4) = 6.12323e-17 (expected: 0)

    === CASE 4: Wrapper for external class ===
    f5(0) = 1 (expected: 1)
    f5(1) = 0.367879 (expected: 0.367879)
*/
}

void Readme_Interpolation()
{
    std::cout << "\n***********************************************************************" << std::endl;
    std::cout << "****                    README - Interpolation                     ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    // Create data points
    Vector<Real> x_data{0, 2.5, 5.0, 7.5, 10.0};
    Vector<Real> y_data{0, 1.5, 4.0, 9.5, 16.0};

    std::cout << "\nData points: ";
    for (int i = 0; i < 5; i++)
        std::cout << "(" << x_data[i] << ", " << y_data[i] << ") ";
    std::cout << std::endl;

    // Create interpolations
    LinearInterpRealFunc   linear_interp(x_data, y_data);
    SplineInterpRealFunc   spline_interp(x_data, y_data);
    PolynomInterpRealFunc  poly_interp(x_data, y_data, 3);   // degree 3

    // Evaluate at test point
    Real x = 3.7;
    std::cout << "\nInterpolation at x = " << x << ":" << std::endl;
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "  Linear:       " << std::setw(12) << linear_interp(x) << std::endl;
    std::cout << "  Cubic Spline: " << std::setw(12) << spline_interp(x) << std::endl;
    std::cout << "  Polynomial:   " << std::setw(12) << poly_interp(x) << std::endl;

    // Compare interpolations at multiple points
    std::cout << "\nComparison across range:" << std::endl;
    std::cout << "    x       Linear    Spline    Polynom" << std::endl;
    std::cout << "  -----    -------   -------   -------" << std::endl;
    for (Real xi = 1.0; xi <= 9.0; xi += 2.0) {
        std::cout << std::setw(6) << xi 
                  << std::setw(11) << linear_interp(xi)
                  << std::setw(10) << spline_interp(xi)
                  << std::setw(10) << poly_interp(xi) << std::endl;
    }

/* Expected OUTPUT:
    Data points: (0, 0) (2.5, 1.5) (5, 4) (7.5, 9.5) (10, 16)

    Interpolation at x = 3.7:
      Linear:          2.700000
      Cubic Spline:    2.409038
      Polynomial:      2.429920

    Comparison across range:
        x       Linear    Spline    Polynom
      -----    -------   -------   -------
       1.0    0.600000  0.464640  0.576000
       3.0    1.900000  1.798080  1.824000
       5.0    4.000000  4.000000  4.000000
       7.0    7.800000  7.742400  7.936000
       9.0   13.400000 13.598400 13.824000
*/
}

void Readme_FunctionsInterpolation()
{
    Readme_DefiningFunctions();
    Readme_Interpolation();
}
