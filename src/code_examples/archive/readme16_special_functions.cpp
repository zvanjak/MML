///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        readme16_special_functions.cpp                                      ///
///  Description: README example - Special Mathematical Functions                     ///
///               Demonstrates Bessel, Legendre, Gamma, Error functions               ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///////////////////////////////////////////////////////////////////////////////////////////

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/StandardFunctions.h"
#endif

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace MML;

void Readme_SpecialFunctions()
{
    std::cout << std::endl;
    std::cout << "=== Special Mathematical Functions ===" << std::endl;
    std::cout << std::fixed << std::setprecision(10);

#if MML_HAS_STD_SPECIAL_FUNCTIONS
    // Bessel Functions (C++17 feature)
    std::cout << std::endl << "--- Bessel Functions ---" << std::endl;
    std::cout << "Bessel functions of the first kind Jn(x):" << std::endl;

    std::cout << std::setw(8) << "x" 
              << std::setw(16) << "J0(x)" 
              << std::setw(16) << "J1(x)" 
              << std::setw(16) << "J2(x)" << std::endl;
    std::cout << std::string(56, '-') << std::endl;

    for (Real x = 0.0; x <= 10.0; x += 2.0)
    {
        std::cout << std::setw(8) << std::setprecision(1) << x
                  << std::setw(16) << std::setprecision(8) << StdFunctions::CylBesselJ(0, x)
                  << std::setw(16) << StdFunctions::CylBesselJ(1, x)
                  << std::setw(16) << StdFunctions::CylBesselJ(2, x) << std::endl;
    }

    // Modified Bessel functions
    std::cout << std::endl << "Modified Bessel functions I0(x), K0(x):" << std::endl;
    std::cout << std::setw(8) << "x" 
              << std::setw(16) << "I0(x)" 
              << std::setw(16) << "K0(x)" << std::endl;
    std::cout << std::string(40, '-') << std::endl;

    for (Real x = 0.5; x <= 3.0; x += 0.5)
    {
        std::cout << std::setw(8) << std::setprecision(1) << x
                  << std::setw(16) << std::setprecision(8) << StdFunctions::CylBesselI(0, x)
                  << std::setw(16) << StdFunctions::CylBesselK(0, x) << std::endl;
    }

    // Legendre Polynomials
    std::cout << std::endl << "--- Legendre Polynomials ---" << std::endl;
    std::cout << "Pn(x) for x = 0.5:" << std::endl;

    Real x_leg = 0.5;
    for (unsigned int n = 0; n <= 5; n++)
    {
        std::cout << "  P" << n << "(" << x_leg << ") = " 
                  << std::setprecision(8) << StdFunctions::Legendre(n, x_leg) << std::endl;
    }

    // Associated Legendre functions
    std::cout << std::endl << "Associated Legendre functions Plm(x) for x = 0.5:" << std::endl;
    std::cout << "  P2,0 = " << StdFunctions::AssocLegendre(2, 0, x_leg) << std::endl;
    std::cout << "  P2,1 = " << StdFunctions::AssocLegendre(2, 1, x_leg) << std::endl;
    std::cout << "  P2,2 = " << StdFunctions::AssocLegendre(2, 2, x_leg) << std::endl;
    std::cout << "  P3,2 = " << StdFunctions::AssocLegendre(3, 2, x_leg) << std::endl;

    // Beta Function
    std::cout << std::endl << "--- Beta Function ---" << std::endl;
    std::cout << "B(a,b) = Γ(a)Γ(b)/Γ(a+b):" << std::endl;
    std::cout << "  B(2, 3) = " << StdFunctions::Beta(2.0, 3.0) 
              << " (should be 1/12 ≈ 0.0833)" << std::endl;
    std::cout << "  B(0.5, 0.5) = " << StdFunctions::Beta(0.5, 0.5) 
              << " (should be π ≈ " << Constants::PI << ")" << std::endl;
#else
    std::cout << "(Note: Full special functions require C++17 std::special_math)" << std::endl;
    std::cout << "Using fallback implementations for Legendre polynomials..." << std::endl;
    
    std::cout << std::endl << "--- Legendre Polynomials (fallback) ---" << std::endl;
    Real x_leg = 0.5;
    for (unsigned int n = 0; n <= 5; n++)
    {
        std::cout << "  P" << n << "(" << x_leg << ") = " 
                  << std::setprecision(8) << StdFunctions::Legendre(n, x_leg) << std::endl;
    }
#endif

    // Gamma Function (always available via std::tgamma)
    std::cout << std::endl << "--- Gamma Function ---" << std::endl;
    std::cout << "Γ(n) = (n-1)! for positive integers:" << std::endl;

    for (int n = 1; n <= 6; n++)
    {
        std::cout << "  Γ(" << n << ") = " << std::setprecision(1) 
                  << StdFunctions::TGamma((Real)n) << std::endl;
    }

    std::cout << std::endl << "Γ(x) for non-integers:" << std::endl;
    std::cout << "  Γ(0.5) = √π = " << std::setprecision(10) 
              << StdFunctions::TGamma(0.5) 
              << " (exact: " << std::sqrt(Constants::PI) << ")" << std::endl;
    std::cout << "  Γ(1.5) = √π/2 = " << StdFunctions::TGamma(1.5) << std::endl;
    std::cout << "  Γ(2.5) = 3√π/4 = " << StdFunctions::TGamma(2.5) << std::endl;

    // Log Gamma (for large arguments)
    std::cout << std::endl << "Log Gamma for large arguments:" << std::endl;
    std::cout << "  ln Γ(100) = " << StdFunctions::LGamma(100.0) << std::endl;
    std::cout << "  ln Γ(1000) = " << StdFunctions::LGamma(1000.0) << std::endl;

    // Error Function (always available)
    std::cout << std::endl << "--- Error Function ---" << std::endl;
    std::cout << "erf(x) and erfc(x) = 1 - erf(x):" << std::endl;
    std::cout << std::setw(8) << "x" 
              << std::setw(16) << "erf(x)" 
              << std::setw(16) << "erfc(x)" << std::endl;
    std::cout << std::string(40, '-') << std::endl;

    for (Real xx = 0.0; xx <= 3.0; xx += 0.5)
    {
        std::cout << std::setw(8) << std::setprecision(1) << xx
                  << std::setw(16) << std::setprecision(10) << StdFunctions::Erf(xx)
                  << std::setw(16) << StdFunctions::Erfc(xx) << std::endl;
    }

    // Verify erf properties
    std::cout << std::endl << "Error function properties:" << std::endl;
    std::cout << "  erf(0) = " << StdFunctions::Erf(0.0) << " (should be 0)" << std::endl;
    std::cout << "  erf(∞) → " << StdFunctions::Erf(10.0) << " (should be 1)" << std::endl;
    std::cout << "  erf(-x) = -erf(x): erf(-1) = " << StdFunctions::Erf(-1.0) 
              << ", -erf(1) = " << -StdFunctions::Erf(1.0) << std::endl;

    std::cout << std::endl;
}
