///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary Docs Demo                              ///
///  File: docs_demo_chebyshev.cpp                                                    ///
///  Demonstrates: ChebyshevApproximation.h - Chebyshev polynomial approximation     ///
///////////////////////////////////////////////////////////////////////////////////////////

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"
#include "algorithms/ChebyshevApproximation.h"
#endif

#include <iostream>
#include <iomanip>
#include <string>

using namespace MML;

namespace {
    inline void DemoHeader(const std::string& name, const std::string& desc) {
        std::cout << "\n========================================\n";
        std::cout << "  " << name << "\n";
        std::cout << "  " << desc << "\n";
        std::cout << "========================================\n\n";
    }
}

namespace MML::Demos::Algorithms
{
    void Chebyshev_BasicApproximation()
    {
        DemoHeader("Chebyshev_BasicApproximation", "Creating Chebyshev approximations");
        
        std::cout << "ChebyshevApproximation represents f(x) as:\n";
        std::cout << "  f(x) ≈ Σ cⱼ Tⱼ(y), where y maps [a,b] to [-1,1]\n\n";
        
        // Approximate sin(x) on [0, π]
        auto sinFunc = [](Real x) { return std::sin(x); };
        ChebyshevApproximation chebSin(sinFunc, 0, Constants::PI, 20);
        
        std::cout << "Approximating sin(x) on [0, π] with 20 terms:\n\n";
        
        std::cout << "x\t\tsin(x)\t\tCheb(x)\t\tError\n";
        std::cout << std::string(60, '-') << "\n";
        
        for (Real x = 0; x <= Constants::PI; x += Constants::PI / 5) {
            Real exact = std::sin(x);
            Real approx = chebSin(x);
            Real error = std::abs(exact - approx);
            std::cout << std::fixed << std::setprecision(4)
                      << x << "\t\t" << exact << "\t\t" << approx << "\t\t"
                      << std::scientific << error << "\n";
        }
        std::cout << std::fixed;
    }
    
    void Chebyshev_Convergence()
    {
        DemoHeader("Chebyshev_Convergence", "Spectral convergence rate");
        
        std::cout << "Chebyshev approximation converges exponentially for smooth functions.\n\n";
        
        auto expFunc = [](Real x) { return std::exp(x); };
        Real a = -1, b = 1;
        Real testPoint = 0.5;
        Real exact = std::exp(testPoint);
        
        std::cout << "Approximating exp(x) on [-1, 1] at x = 0.5:\n\n";
        std::cout << "n terms\t\tApprox\t\t\tError\n";
        std::cout << std::string(50, '-') << "\n";
        
        for (int n : {5, 10, 15, 20, 30}) {
            ChebyshevApproximation cheb(expFunc, a, b, n);
            Real approx = cheb(testPoint);
            Real error = std::abs(approx - exact);
            std::cout << n << "\t\t" << std::setprecision(12) << approx 
                      << "\t" << std::scientific << error << std::fixed << "\n";
        }
        
        std::cout << "\nExact value: " << std::setprecision(12) << exact << "\n";
    }
    
    void Chebyshev_Truncation()
    {
        DemoHeader("Chebyshev_Truncation", "Truncation for economization");
        
        std::cout << "Use Truncate(m) to use only first m terms, reducing computation.\n\n";
        
        auto cosFunc = [](Real x) { return std::cos(x); };
        ChebyshevApproximation cheb(cosFunc, -Constants::PI, Constants::PI, 30);
        
        std::cout << "cos(x) on [-π, π] with 30 coefficients:\n\n";
        
        Real testPoint = 1.0;
        Real exact = std::cos(testPoint);
        
        std::cout << "Testing at x = 1.0 (exact = " << exact << "):\n\n";
        std::cout << "m terms\t\tApprox\t\t\tError\n";
        std::cout << std::string(50, '-') << "\n";
        
        for (int m : {5, 10, 15, 20, 25, 30}) {
            cheb.Truncate(m);
            Real approx = cheb(testPoint);
            Real error = std::abs(approx - exact);
            std::cout << m << "\t\t" << std::setprecision(10) << approx 
                      << "\t" << std::scientific << error << std::fixed << "\n";
        }
    }
    
    void Chebyshev_Derivative()
    {
        DemoHeader("Chebyshev_Derivative", "Computing derivatives");
        
        std::cout << "Derivative() returns a new ChebyshevApproximation of f'(x).\n\n";
        
        // f(x) = sin(x), f'(x) = cos(x)
        auto sinFunc = [](Real x) { return std::sin(x); };
        ChebyshevApproximation chebSin(sinFunc, 0, Constants::PI, 25);
        
        ChebyshevApproximation chebDerivative = chebSin.Derivative();
        
        std::cout << "f(x) = sin(x), computing f'(x) via Chebyshev:\n\n";
        
        std::cout << "x\t\tcos(x) exact\tf'(x) Cheb\tError\n";
        std::cout << std::string(60, '-') << "\n";
        
        for (Real x = 0.1; x < Constants::PI; x += 0.5) {
            Real exact = std::cos(x);
            Real approx = chebDerivative(x);
            Real error = std::abs(exact - approx);
            std::cout << std::fixed << std::setprecision(4)
                      << x << "\t\t" << exact << "\t\t" << approx << "\t\t"
                      << std::scientific << error << "\n";
        }
        std::cout << std::fixed;
    }
    
    void Chebyshev_Integral()
    {
        DemoHeader("Chebyshev_Integral", "Computing integrals");
        
        std::cout << "Integral() returns a new ChebyshevApproximation of ∫f(x)dx.\n\n";
        
        // f(x) = x², ∫f = x³/3
        auto sqFunc = [](Real x) { return x * x; };
        ChebyshevApproximation chebSq(sqFunc, 0, 2, 20);
        
        ChebyshevApproximation chebIntegral = chebSq.Integral();
        
        std::cout << "f(x) = x², computing ∫f(x)dx via Chebyshev:\n";
        std::cout << "Exact: ∫x² dx = x³/3 (plus constant)\n\n";
        
        std::cout << "x\t\tx³/3 - C\t∫f Cheb\t\tDifference\n";
        std::cout << std::string(60, '-') << "\n";
        
        // The integral has an arbitrary constant; compare differences
        Real c = chebIntegral(0) - 0.0;  // Constant offset
        
        for (Real x = 0; x <= 2; x += 0.4) {
            Real exact = x * x * x / 3.0;
            Real approx = chebIntegral(x) - c;
            Real diff = std::abs(exact - approx);
            std::cout << std::fixed << std::setprecision(4)
                      << x << "\t\t" << exact << "\t\t" << approx << "\t\t"
                      << std::scientific << diff << "\n";
        }
        std::cout << std::fixed;
    }
    
    void Chebyshev_RootFinding()
    {
        DemoHeader("Chebyshev_RootFinding", "Chebyshev polynomial roots");
        
        std::cout << "ChebyshevRoots(n) returns the n roots of Chebyshev polynomial T_n.\n";
        std::cout << "These roots are used as interpolation nodes for optimal approximation.\n\n";
        
        // Get Chebyshev nodes (roots of T_n)
        int n = 5;
        Vector<Real> roots = ChebyshevRoots(n);
        
        std::cout << "Roots of T_" << n << " (Chebyshev nodes):\n";
        for (int i = 0; i < n; ++i) {
            std::cout << "  x_" << i << " = " << std::fixed << std::setprecision(6) 
                      << roots[i] << "\n";
        }
        
        std::cout << "\nThese are cos((2k+1)π/(2n)) for k = 0, 1, ..., n-1\n";
        std::cout << "They minimize interpolation error (Runge phenomenon).\n";
    }
    
    void Chebyshev_Applications()
    {
        DemoHeader("Chebyshev_Applications", "Practical applications");
        
        std::cout << "Common uses for Chebyshev approximation:\n\n";
        
        std::cout << "1. FUNCTION TABULATION\n";
        std::cout << "   Replace expensive function evaluation with polynomial.\n";
        std::cout << "   Example: Precompute Chebyshev coeffs for special functions.\n\n";
        
        std::cout << "2. NUMERICAL INTEGRATION\n";
        std::cout << "   Clenshaw-Curtis quadrature uses Chebyshev nodes.\n";
        std::cout << "   Achieves spectral convergence for smooth integrands.\n\n";
        
        std::cout << "3. DIFFERENTIAL EQUATIONS\n";
        std::cout << "   Spectral methods represent solutions as Chebyshev series.\n";
        std::cout << "   High accuracy with relatively few terms.\n\n";
        
        std::cout << "4. FUNCTION ECONOMIZATION\n";
        std::cout << "   Truncate series to desired accuracy.\n";
        std::cout << "   Chebyshev is near-minimax optimal.\n\n";
        
        std::cout << "5. ROOT FINDING\n";
        std::cout << "   Convert to Chebyshev, find roots via eigenvalue problem.\n";
        std::cout << "   Robust and accurate.\n";
    }
}

void Docs_Demo_ChebyshevApproximation()
{
    using namespace MML::Demos::Algorithms;
    
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "               CHEBYSHEV APPROXIMATION DEMO\n";
    std::cout << "               ChebyshevApproximation.h\n";
    std::cout << std::string(70, '=') << "\n";
    std::cout << "This demo covers:\n";
    std::cout << "  - Creating Chebyshev approximations from functions\n";
    std::cout << "  - Spectral convergence for smooth functions\n";
    std::cout << "  - Truncation for economization\n";
    std::cout << "  - Derivative computation\n";
    std::cout << "  - Integration\n";
    std::cout << "  - Root finding\n";
    std::cout << std::string(70, '=') << "\n\n";
    
    Chebyshev_BasicApproximation();
    Chebyshev_Convergence();
    Chebyshev_Truncation();
    Chebyshev_Derivative();
    Chebyshev_Integral();
    Chebyshev_RootFinding();
    Chebyshev_Applications();
}
