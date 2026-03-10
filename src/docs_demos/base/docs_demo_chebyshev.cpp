///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        docs_demo_chebyshev.cpp                                             ///
///  Description: Documentation examples for Chebyshev polynomials                    ///
///               These examples are referenced from docs/base/ChebyshevPolynom.md    ///
///                                                                                   ///
///  Usage:       Run MML_DocsDemo application to execute these examples              ///
///////////////////////////////////////////////////////////////////////////////////////////

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/ChebyshevPolynom.h"
#endif

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace MML;

///////////////////////////////////////////////////////////////////////////////////////////
///                      Chebyshev Polynomials of the First Kind                        ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_Chebyshev_FirstKind()
{
    std::cout << "\n=== Chebyshev Polynomials of the First Kind T_n(x) ===\n\n";

    // Evaluate T_n at specific points
    std::cout << "T_n(x) evaluation table:" << std::endl;
    std::cout << std::setw(8) << "x" << " | ";
    for (int n = 0; n <= 5; ++n) {
        std::cout << std::setw(10) << "T_" << n << "(x)";
    }
    std::cout << std::endl;
    std::cout << std::string(75, '-') << std::endl;

    double xValues[] = {-1.0, -0.5, 0.0, 0.5, 1.0};
    for (double x : xValues) {
        std::cout << std::fixed << std::setprecision(2) << std::setw(8) << x << " | ";
        for (int n = 0; n <= 5; ++n) {
            std::cout << std::setw(12) << std::setprecision(4) << ChebyshevT(n, x);
        }
        std::cout << std::endl;
    }

    // Key property: T_n(1) = 1 for all n
    std::cout << "\nProperty: T_n(1) = 1 for all n" << std::endl;
    std::cout << "  ";
    for (int n = 0; n <= 6; ++n) {
        std::cout << "T_" << n << "(1)=" << ChebyshevT(n, 1.0) << "  ";
    }
    std::cout << std::endl;

    // Key property: T_n(-1) = (-1)^n
    std::cout << "\nProperty: T_n(-1) = (-1)^n" << std::endl;
    std::cout << "  ";
    for (int n = 0; n <= 6; ++n) {
        std::cout << "T_" << n << "(-1)=" << ChebyshevT(n, -1.0) << "  ";
    }
    std::cout << std::endl;

    // Trigonometric relationship: T_n(cos θ) = cos(nθ)
    std::cout << "\nTrigonometric identity: T_n(cos θ) = cos(nθ)" << std::endl;
    double theta = Constants::PI / 6.0;  // 30 degrees
    double cosTheta = std::cos(theta);
    std::cout << "  θ = π/6 (30°), cos(θ) = " << std::fixed << std::setprecision(4) << cosTheta << std::endl;
    
    for (int n = 0; n <= 4; ++n) {
        double chebyshev = ChebyshevT(n, cosTheta);
        double cosnTheta = std::cos(n * theta);
        std::cout << "  T_" << n << "(cos θ) = " << std::setw(8) << chebyshev 
                  << ",  cos(" << n << "θ) = " << std::setw(8) << cosnTheta
                  << ",  diff = " << std::scientific << std::setprecision(2) 
                  << std::abs(chebyshev - cosnTheta) << std::endl;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
///                      Chebyshev Polynomials of the Second Kind                       ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_Chebyshev_SecondKind()
{
    std::cout << "\n=== Chebyshev Polynomials of the Second Kind U_n(x) ===\n\n";

    // Evaluate U_n at specific points
    std::cout << "U_n(x) evaluation table:" << std::endl;
    std::cout << std::setw(8) << "x" << " | ";
    for (int n = 0; n <= 5; ++n) {
        std::cout << std::setw(10) << "U_" << n << "(x)";
    }
    std::cout << std::endl;
    std::cout << std::string(75, '-') << std::endl;

    double xValues[] = {-1.0, -0.5, 0.0, 0.5, 1.0};
    for (double x : xValues) {
        std::cout << std::fixed << std::setprecision(2) << std::setw(8) << x << " | ";
        for (int n = 0; n <= 5; ++n) {
            std::cout << std::setw(12) << std::setprecision(4) << ChebyshevU(n, x);
        }
        std::cout << std::endl;
    }

    // Key property: U_n(1) = n + 1
    std::cout << "\nProperty: U_n(1) = n + 1" << std::endl;
    std::cout << "  ";
    for (int n = 0; n <= 6; ++n) {
        std::cout << "U_" << n << "(1)=" << ChebyshevU(n, 1.0) << "  ";
    }
    std::cout << std::endl;

    // Key property: U_n(-1) = (-1)^n * (n + 1)
    std::cout << "\nProperty: U_n(-1) = (-1)^n * (n + 1)" << std::endl;
    std::cout << "  ";
    for (int n = 0; n <= 6; ++n) {
        std::cout << "U_" << n << "(-1)=" << ChebyshevU(n, -1.0) << "  ";
    }
    std::cout << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
///                           Comparison of T_n and U_n                                 ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_Chebyshev_Comparison()
{
    std::cout << "\n=== Comparison of T_n(x) and U_n(x) ===\n\n";

    std::cout << "Explicit formulas for low degrees:" << std::endl;
    std::cout << "  T_0(x) = 1" << std::endl;
    std::cout << "  T_1(x) = x" << std::endl;
    std::cout << "  T_2(x) = 2x² - 1" << std::endl;
    std::cout << "  T_3(x) = 4x³ - 3x" << std::endl;
    std::cout << "  T_4(x) = 8x⁴ - 8x² + 1" << std::endl;
    std::cout << std::endl;
    std::cout << "  U_0(x) = 1" << std::endl;
    std::cout << "  U_1(x) = 2x" << std::endl;
    std::cout << "  U_2(x) = 4x² - 1" << std::endl;
    std::cout << "  U_3(x) = 8x³ - 4x" << std::endl;
    std::cout << "  U_4(x) = 16x⁴ - 12x² + 1" << std::endl;
    std::cout << std::endl;

    // Verify formulas at x = 0.5
    double x = 0.5;
    std::cout << "Verification at x = 0.5:" << std::endl;
    std::cout << "  T_2(0.5) = 2(0.5)² - 1 = " << (2*x*x - 1) 
              << ", ChebyshevT(2, 0.5) = " << ChebyshevT(2, x) << std::endl;
    std::cout << "  T_3(0.5) = 4(0.5)³ - 3(0.5) = " << (4*x*x*x - 3*x)
              << ", ChebyshevT(3, 0.5) = " << ChebyshevT(3, x) << std::endl;
    std::cout << "  U_2(0.5) = 4(0.5)² - 1 = " << (4*x*x - 1)
              << ", ChebyshevU(2, 0.5) = " << ChebyshevU(2, x) << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
///                           Orthogonality Property                                    ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_Chebyshev_Orthogonality()
{
    std::cout << "\n=== Orthogonality Properties ===\n\n";

    std::cout << "Chebyshev polynomials are orthogonal with respect to different weights:" << std::endl;
    std::cout << std::endl;
    std::cout << "T_n orthogonality (weight w(x) = 1/√(1-x²)):" << std::endl;
    std::cout << "  ∫₋₁¹ T_m(x) T_n(x) / √(1-x²) dx = {0 if m≠n, π if m=n=0, π/2 if m=n>0}" << std::endl;
    std::cout << std::endl;
    std::cout << "U_n orthogonality (weight w(x) = √(1-x²)):" << std::endl;
    std::cout << "  ∫₋₁¹ U_m(x) U_n(x) √(1-x²) dx = {0 if m≠n, π/2 if m=n}" << std::endl;
    std::cout << std::endl;

    // Demonstrate zeros of T_n(x)
    std::cout << "Zeros of T_n(x):" << std::endl;
    std::cout << "  The n zeros of T_n are: x_k = cos((2k-1)π/(2n)), k = 1,...,n" << std::endl;
    std::cout << std::endl;
    
    std::cout << "  Zeros of T_4:" << std::endl;
    for (int k = 1; k <= 4; ++k) {
        double x_k = std::cos((2.0*k - 1.0) * Constants::PI / 8.0);
        double T4_at_xk = ChebyshevT(4, x_k);
        std::cout << "    x_" << k << " = cos(" << (2*k-1) << "π/8) = " 
                  << std::fixed << std::setprecision(6) << x_k
                  << ",  T_4(x_" << k << ") = " << std::scientific << std::setprecision(2) 
                  << T4_at_xk << std::endl;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
///                           Applications                                              ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_Chebyshev_Applications()
{
    std::cout << "\n=== Applications of Chebyshev Polynomials ===\n\n";

    std::cout << "1. FUNCTION APPROXIMATION (Chebyshev expansion)" << std::endl;
    std::cout << "   Any continuous function on [-1,1] can be approximated as:" << std::endl;
    std::cout << "   f(x) ≈ c_0/2 + Σ c_n T_n(x)" << std::endl;
    std::cout << "   with coefficients c_n = (2/π) ∫₋₁¹ f(x)T_n(x)/√(1-x²) dx" << std::endl;
    std::cout << std::endl;

    std::cout << "2. CLENSHAW-CURTIS QUADRATURE" << std::endl;
    std::cout << "   Integration using Chebyshev nodes (zeros of T_n)" << std::endl;
    std::cout << "   Nodes: x_k = cos((2k-1)π/(2n)), k = 1,...,n" << std::endl;
    std::cout << std::endl;

    std::cout << "3. SPECTRAL METHODS FOR PDEs" << std::endl;
    std::cout << "   Chebyshev collocation for solving differential equations" << std::endl;
    std::cout << "   Provides spectral accuracy for smooth solutions" << std::endl;
    std::cout << std::endl;

    std::cout << "4. MINIMAX APPROXIMATION" << std::endl;
    std::cout << "   T_n(x)/2^(n-1) is the monic polynomial with smallest max error on [-1,1]" << std::endl;
    std::cout << "   This is the equioscillation property" << std::endl;
    std::cout << std::endl;

    // Demonstrate equioscillation
    std::cout << "Demonstrating equioscillation for T_4 / 8 (monic degree-4 polynomial):" << std::endl;
    double extrema[] = {-1.0, -0.7071, 0.0, 0.7071, 1.0};  // cos(kπ/4)
    for (double x : extrema) {
        double val = ChebyshevT(4, x) / 8.0;
        std::cout << "  At x = " << std::setw(7) << std::fixed << std::setprecision(4) << x 
                  << ": T_4/8 = " << std::setw(10) << val << std::endl;
    }
    std::cout << "  (Notice: alternates between +1/8 and -1/8 = ±0.125)" << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
///                           Main Entry Point                                          ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_Chebyshev()
{
    std::cout << "╔════════════════════════════════════════════════════════════════════╗" << std::endl;
    std::cout << "║          Chebyshev Polynomial Documentation Examples              ║" << std::endl;
    std::cout << "║               (from docs/base/ChebyshevPolynom.md)                ║" << std::endl;
    std::cout << "╚════════════════════════════════════════════════════════════════════╝" << std::endl;

    Demo_Chebyshev_FirstKind();
    Demo_Chebyshev_SecondKind();
    Demo_Chebyshev_Comparison();
    Demo_Chebyshev_Orthogonality();
    Demo_Chebyshev_Applications();

    std::cout << "\n=== All Chebyshev Examples Complete ===\n" << std::endl;
}
