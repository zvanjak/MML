///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        readme09_polynomials.cpp                                            ///
///  Description: README example - Polynomials & Algebra                              ///
///               Demonstrates polynomial creation, arithmetic, and root finding      ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///////////////////////////////////////////////////////////////////////////////////////////

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/Polynom.h"
#include "algorithms/RootFinding.h"
#endif

#include <iostream>
#include <iomanip>

using namespace MML;

void Readme_Polynomials()
{
    std::cout << std::endl;
    std::cout << "=== Polynomials & Algebra ===" << std::endl;

    // Create polynomials: p(x) = 2x³ - 3x² + x - 5
    PolynomReal p{-5, 1, -3, 2};  // Coefficients: [constant, x, x², x³]

    // Evaluate at a point
    Real val = p(2.0);  // p(2) = 2(8) - 3(4) + 2 - 5 = 16 - 12 + 2 - 5 = 1
    std::cout << "p(x) = 2x³ - 3x² + x - 5" << std::endl;
    std::cout << "p(2) = " << val << std::endl;

    // Polynomial arithmetic
    PolynomReal q{1, 2};           // q(x) = 2x + 1
    PolynomReal sum = p + q;       // Addition
    PolynomReal prod = p * q;      // Multiplication
    std::cout << "q(x) = 2x + 1" << std::endl;
    std::cout << "Degree of p*q: " << prod.degree() << std::endl;

    // Calculus on polynomials
    PolynomReal dp = p.derivative();        // p'(x) = 6x² - 6x + 1
    PolynomReal ip = p.integral();          // ∫p(x)dx (constant = 0)
    std::cout << "p'(x) at x=1: " << dp(1.0) << std::endl;

    // Solve quadratic: x² - 5x + 6 = 0 → x = 2, 3
    Complex r1, r2;
    int numReal = SolveQuadratic(1.0, -5.0, 6.0, r1, r2);
    std::cout << std::endl << "Quadratic x² - 5x + 6 = 0:" << std::endl;
    std::cout << "  Roots: " << r1.real() << ", " << r2.real() << std::endl;

    // Solve cubic: x³ - 6x² + 11x - 6 = 0 → x = 1, 2, 3
    Complex c1, c2, c3;
    SolveCubic(1.0, -6.0, 11.0, -6.0, c1, c2, c3);
    std::cout << "Cubic roots: " << c1.real() << ", " << c2.real() 
              << ", " << c3.real() << std::endl;

    std::cout << std::endl;
}
