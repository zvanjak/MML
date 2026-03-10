///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        docs_demo_curve_fitting.cpp                                         ///
///  Description: Demonstration of CurveFitting.h functionality                       ///
///               - Linear least squares fitting                                      ///
///               - Polynomial fitting                                                ///
///               - General linear least squares with custom basis functions          ///
///               - Fitting statistics and quality metrics                            ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                    ///
///////////////////////////////////////////////////////////////////////////////////////////

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"
#include "base/Vector/Vector.h"
#include "algorithms/CurveFitting.h"
#endif

#include <iostream>
#include <iomanip>
#include <cmath>
#include <random>

using namespace MML;

///////////////////////////////////////////////////////////////////////////////////////////
///                              LINEAR LEAST SQUARES                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_LinearLeastSquares()
{
    std::cout << "\n=== LINEAR LEAST SQUARES ===\n" << std::endl;
    
    // Generate data points along y = 2x + 1 with some noise
    Vector<Real> x_data({0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0});
    Vector<Real> y_data({1.2, 2.9, 5.1, 7.0, 8.8, 11.2, 12.9, 15.1, 17.0, 18.8});
    
    std::cout << "Data points (noisy line y ≈ 2x + 1):" << std::endl;
    std::cout << "  x: "; x_data.Print(std::cout, 5, 1);
    std::cout << "  y: "; y_data.Print(std::cout, 5, 1);
    
    // Simple linear fit
    Real a, b;
    Real residual = LinearLeastSquares(x_data, y_data, a, b);
    
    std::cout << "\nLinear fit: y = " << a << " * x + " << b << std::endl;
    std::cout << "Residual norm: " << residual << std::endl;
    std::cout << "Expected: y ≈ 2x + 1" << std::endl;
    
    // Detailed linear fit with statistics
    std::cout << "\n--- Detailed Linear Fit ---" << std::endl;
    LinearFitResult<Real> detailed = LinearLeastSquaresDetailed(x_data, y_data);
    
    std::cout << "Slope (a): " << detailed.a << std::endl;
    std::cout << "Intercept (b): " << detailed.b << std::endl;
    std::cout << "R² (coefficient of determination): " << detailed.r_squared << std::endl;
    std::cout << "Mean squared error: " << detailed.mean_squared_error << std::endl;
    std::cout << "\n(R² = 1.0 means perfect fit, R² = 0.0 means no linear relationship)" << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
///                             POLYNOMIAL FITTING                                      ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_PolynomialFit()
{
    std::cout << "\n=== POLYNOMIAL FITTING ===\n" << std::endl;
    
    // Generate data from y = x² - 2x + 1 = (x-1)² with noise
    Vector<Real> x_data({-2.0, -1.0, 0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0});
    Vector<Real> y_data({9.1,   4.0, 1.2, 0.3, 0.1, 0.2, 1.1, 4.0, 9.0}); // y ≈ (x-1)²
    
    std::cout << "Data points (noisy parabola y ≈ x² - 2x + 1):" << std::endl;
    std::cout << "  x: "; x_data.Print(std::cout, 6, 1);
    std::cout << "  y: "; y_data.Print(std::cout, 6, 1);
    
    // Fit polynomials of different degrees
    for (int degree = 1; degree <= 4; degree++)
    {
        GeneralLinearFitResult<Real> result = PolynomialFit(x_data, y_data, degree);
        
        std::cout << "\n--- Degree " << degree << " Polynomial ---" << std::endl;
        std::cout << "Coefficients: ";
        for (int i = 0; i < result.coefficients.size(); i++) {
            if (i > 0 && result.coefficients[i] >= 0) std::cout << "+ ";
            std::cout << result.coefficients[i];
            if (i == 1) std::cout << "*x ";
            else if (i > 1) std::cout << "*x^" << i << " ";
        }
        std::cout << std::endl;
        std::cout << "R²: " << result.r_squared << std::endl;
        std::cout << "Adjusted R²: " << result.adjusted_r_squared << std::endl;
        std::cout << "Residual norm: " << result.residual_norm << std::endl;
    }
    
    std::cout << "\n(Degree 2 should give best fit with R² ≈ 1)" << std::endl;
    std::cout << "(Higher degrees may overfit with similar or worse R²)" << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
///                     EVALUATING POLYNOMIAL FIT                                       ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_EvaluatePolynomial()
{
    std::cout << "\n=== EVALUATING POLYNOMIAL FIT ===\n" << std::endl;
    
    // Fit a cubic polynomial
    Vector<Real> x_data({0.0, 1.0, 2.0, 3.0, 4.0, 5.0});
    Vector<Real> y_data({1.0, 2.0, 5.0, 14.0, 35.0, 74.0}); // y = x³ - x + 1 + noise
    
    GeneralLinearFitResult<Real> fit = PolynomialFit(x_data, y_data, 3);
    
    std::cout << "Fitted cubic polynomial coefficients:" << std::endl;
    std::cout << "  c₀ = " << fit.coefficients[0] << std::endl;
    std::cout << "  c₁ = " << fit.coefficients[1] << std::endl;
    std::cout << "  c₂ = " << fit.coefficients[2] << std::endl;
    std::cout << "  c₃ = " << fit.coefficients[3] << std::endl;
    
    std::cout << "\nEvaluating fitted polynomial at various points:" << std::endl;
    std::cout << std::setw(8) << "x" << std::setw(15) << "y_data" 
              << std::setw(15) << "y_fitted" << std::setw(12) << "error" << std::endl;
    std::cout << std::string(50, '-') << std::endl;
    
    for (int i = 0; i < x_data.size(); i++) {
        Real x = x_data[i];
        Real y_fitted = EvaluatePolynomial(x, fit.coefficients);
        Real error = y_data[i] - y_fitted;
        std::cout << std::fixed << std::setprecision(2)
                  << std::setw(8) << x 
                  << std::setw(15) << y_data[i]
                  << std::setw(15) << y_fitted
                  << std::setw(12) << error << std::endl;
    }
    
    // Interpolate at points not in the data
    std::cout << "\nInterpolating at new points:" << std::endl;
    for (Real x : {0.5, 1.5, 2.5, 3.5, 4.5}) {
        Real y = EvaluatePolynomial(x, fit.coefficients);
        std::cout << "  f(" << x << ") = " << y << std::endl;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
///                    GENERAL LINEAR LEAST SQUARES                                     ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_GeneralLinearLeastSquares()
{
    std::cout << "\n=== GENERAL LINEAR LEAST SQUARES ===\n" << std::endl;
    
    std::cout << "Fitting y = c₀ + c₁*sin(x) + c₂*cos(x) + c₃*sin(2x) to data\n" << std::endl;
    
    // Generate data from y = 1 + 2*sin(x) + 0.5*cos(x)
    Vector<Real> x_data(20);
    Vector<Real> y_data(20);
    
    // Use a random number generator for noise
    std::mt19937 gen(42);  // Fixed seed for reproducibility
    std::normal_distribution<Real> noise(0.0, 0.1);
    
    for (int i = 0; i < 20; i++) {
        x_data[i] = i * 0.5;
        y_data[i] = 1.0 + 2.0*std::sin(x_data[i]) + 0.5*std::cos(x_data[i]) + noise(gen);
    }
    
    std::cout << "First 10 data points:" << std::endl;
    std::cout << "  x: ";
    for (int i = 0; i < 10; i++) std::cout << std::fixed << std::setprecision(2) << x_data[i] << " ";
    std::cout << std::endl;
    std::cout << "  y: ";
    for (int i = 0; i < 10; i++) std::cout << std::fixed << std::setprecision(2) << y_data[i] << " ";
    std::cout << std::endl;
    
    // Define basis functions using lambdas
    Vector<std::function<Real(Real)>> basis(4);
    basis[0] = [](Real) -> Real { return 1.0; };           // Constant
    basis[1] = [](Real x) -> Real { return std::sin(x); }; // sin(x)
    basis[2] = [](Real x) -> Real { return std::cos(x); }; // cos(x)
    basis[3] = [](Real x) -> Real { return std::sin(2*x); }; // sin(2x)
    
    GeneralLinearFitResult<Real> result = GeneralLinearLeastSquares(x_data, y_data, basis);
    
    std::cout << "\nFitted coefficients:" << std::endl;
    std::cout << "  c₀ (constant): " << result.coefficients[0] << " (expected: 1.0)" << std::endl;
    std::cout << "  c₁ (sin(x)):   " << result.coefficients[1] << " (expected: 2.0)" << std::endl;
    std::cout << "  c₂ (cos(x)):   " << result.coefficients[2] << " (expected: 0.5)" << std::endl;
    std::cout << "  c₃ (sin(2x)):  " << result.coefficients[3] << " (expected: 0.0)" << std::endl;
    
    std::cout << "\nFit statistics:" << std::endl;
    std::cout << "  R²: " << result.r_squared << std::endl;
    std::cout << "  Adjusted R²: " << result.adjusted_r_squared << std::endl;
    std::cout << "  Residual norm: " << result.residual_norm << std::endl;
    std::cout << "  Mean squared error: " << result.mean_squared_error << std::endl;
    std::cout << "  Effective rank: " << result.effective_rank << std::endl;
    std::cout << "  Condition number: " << result.condition_number << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        EXPONENTIAL FITTING                                          ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_ExponentialFit()
{
    std::cout << "\n=== EXPONENTIAL FITTING ===\n" << std::endl;
    
    std::cout << "Fitting y = c₀ + c₁*exp(-x) + c₂*exp(-2x) to decay data\n" << std::endl;
    
    // Generate data from y = 1 + 3*exp(-x) + 0.5*exp(-2x)
    Vector<Real> x_data({0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0});
    Vector<Real> y_data(x_data.size());
    
    for (int i = 0; i < x_data.size(); i++) {
        Real x = x_data[i];
        y_data[i] = 1.0 + 3.0*std::exp(-x) + 0.5*std::exp(-2*x);
    }
    
    std::cout << "Data (exact, no noise):" << std::endl;
    std::cout << "  x: "; x_data.Print(std::cout, 5, 1);
    std::cout << "  y: "; y_data.Print(std::cout, 6, 3);
    
    // Define exponential basis functions
    Vector<std::function<Real(Real)>> basis(3);
    basis[0] = [](Real) -> Real { return 1.0; };
    basis[1] = [](Real x) -> Real { return std::exp(-x); };
    basis[2] = [](Real x) -> Real { return std::exp(-2*x); };
    
    GeneralLinearFitResult<Real> result = GeneralLinearLeastSquares(x_data, y_data, basis);
    
    std::cout << "\nFitted coefficients:" << std::endl;
    std::cout << "  c₀ (constant):  " << result.coefficients[0] << " (expected: 1.0)" << std::endl;
    std::cout << "  c₁ (exp(-x)):   " << result.coefficients[1] << " (expected: 3.0)" << std::endl;
    std::cout << "  c₂ (exp(-2x)):  " << result.coefficients[2] << " (expected: 0.5)" << std::endl;
    std::cout << "\nR²: " << result.r_squared << " (should be ≈ 1.0 for exact data)" << std::endl;
    
    // Evaluate the fit
    std::cout << "\nEvaluating fit at x = 2.5:" << std::endl;
    Real x_test = 2.5;
    Real y_exact = 1.0 + 3.0*std::exp(-x_test) + 0.5*std::exp(-2*x_test);
    Real y_fit = result.evaluate(x_test, basis);
    std::cout << "  Exact: " << y_exact << std::endl;
    std::cout << "  Fitted: " << y_fit << std::endl;
    std::cout << "  Error: " << std::abs(y_exact - y_fit) << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
///                          EDGE CASES AND ROBUSTNESS                                  ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_EdgeCases()
{
    std::cout << "\n=== EDGE CASES AND ROBUSTNESS ===\n" << std::endl;
    
    // Two points - exact fit
    std::cout << "--- Two Points (Exact Fit) ---" << std::endl;
    Vector<Real> x2({1.0, 3.0});
    Vector<Real> y2({2.0, 6.0});
    Real a, b;
    Real residual = LinearLeastSquares(x2, y2, a, b);
    std::cout << "Points: (1, 2) and (3, 6)" << std::endl;
    std::cout << "Fit: y = " << a << "*x + " << b << std::endl;
    std::cout << "Residual: " << residual << " (should be 0 for 2 points)" << std::endl;
    
    // Single point
    std::cout << "\n--- Single Point ---" << std::endl;
    Vector<Real> x1({5.0});
    Vector<Real> y1({10.0});
    residual = LinearLeastSquares(x1, y1, a, b);
    std::cout << "Point: (5, 10)" << std::endl;
    std::cout << "Fit: y = " << a << "*x + " << b << " (horizontal line through point)" << std::endl;
    
    // Ill-conditioned polynomial fit (Runge's phenomenon)
    std::cout << "\n--- High-Degree Polynomial Warning ---" << std::endl;
    std::cout << "High-degree polynomials can have large condition numbers" << std::endl;
    std::cout << "and suffer from Runge's phenomenon (oscillation at edges)." << std::endl;
    
    Vector<Real> x_runge({-1.0, -0.5, 0.0, 0.5, 1.0});
    Vector<Real> y_runge({0.5, 0.8, 1.0, 0.8, 0.5}); // Bell-shaped
    
    for (int degree = 2; degree <= 4; degree++) {
        GeneralLinearFitResult<Real> result = PolynomialFit(x_runge, y_runge, degree);
        std::cout << "  Degree " << degree << ": condition number = " 
                  << std::scientific << result.condition_number << std::fixed << std::endl;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
///                               MAIN ENTRY POINT                                      ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_CurveFitting()
{
    std::cout << "###################################################################" << std::endl;
    std::cout << "###               CurveFitting.h - DEMONSTRATION                ###" << std::endl;
    std::cout << "###################################################################" << std::endl;
    
    Demo_LinearLeastSquares();
    Demo_PolynomialFit();
    Demo_EvaluatePolynomial();
    Demo_GeneralLinearLeastSquares();
    Demo_ExponentialFit();
    Demo_EdgeCases();
    
    std::cout << "\n###################################################################" << std::endl;
    std::cout << "###                   DEMONSTRATION COMPLETE                    ###" << std::endl;
    std::cout << "###################################################################\n" << std::endl;
}
