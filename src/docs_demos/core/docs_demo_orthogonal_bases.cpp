///////////////////////////////////////////////////////////////////////////////////////////
// docs_demo_orthogonal_bases.cpp - Demonstrations for OrthogonalBases.md
// 
// This file contains runnable examples for all orthogonal polynomial bases:
//   - Legendre polynomials
//   - Hermite polynomials  
//   - Laguerre polynomials
//   - Chebyshev polynomials (first kind)
//
// Each demo verifies that the documented API is correct and working.
///////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <iomanip>
#include <cmath>
#include <functional>

#include "MMLBase.h"
#include "base/Function.h"
#include "core/OrthogonalBasis/LegendreBasis.h"
#include "core/OrthogonalBasis/HermiteBasis.h"
#include "core/OrthogonalBasis/LaguerreBasis.h"
#include "core/OrthogonalBasis/ChebyshevBasis.h"

using namespace MML;

///////////////////////////////////////////////////////////////////////////////////////////
// Demo: Legendre Polynomials
// Domain: [-1, 1], Weight: w(x) = 1
///////////////////////////////////////////////////////////////////////////////////////////
void Docs_Demo_LegendreBasis()
{
    std::cout << "========================================\n";
    std::cout << "   LEGENDRE POLYNOMIALS DEMO\n";
    std::cout << "========================================\n\n";

    LegendreBasis legendre;

    // Basic evaluation - first few Legendre polynomials
    std::cout << "--- Polynomial Values at x = 0.5 ---\n";
    Real x = 0.5;
    for (int n = 0; n <= 5; ++n)
    {
        std::cout << "P_" << n << "(0.5) = " << std::setw(12) << std::setprecision(8) 
                  << legendre.Evaluate(n, x) << "\n";
    }

    // Known values check
    std::cout << "\n--- Known Value Verification ---\n";
    std::cout << "P_0(x) = 1:       P_0(0.3) = " << legendre.Evaluate(0, 0.3) << "\n";
    std::cout << "P_1(x) = x:       P_1(0.3) = " << legendre.Evaluate(1, 0.3) << "\n";
    std::cout << "P_2(x) = (3x²-1)/2: P_2(0.3) = " << legendre.Evaluate(2, 0.3);
    std::cout << " (expected: " << (3.0*0.3*0.3 - 1.0)/2.0 << ")\n";

    // Domain and weight
    std::cout << "\n--- Domain and Weight ---\n";
    std::cout << "Domain: [" << legendre.DomainMin() << ", " << legendre.DomainMax() << "]\n";
    std::cout << "Weight w(0.5) = " << legendre.WeightFunction(0.5) << " (should be 1.0)\n";

    // Normalization
    std::cout << "\n--- Normalization ||Pₙ||² = 2/(2n+1) ---\n";
    for (int n = 0; n <= 4; ++n)
    {
        Real expected = 2.0 / (2.0 * n + 1.0);
        std::cout << "||P_" << n << "||² = " << std::setw(10) << legendre.Normalization(n)
                  << "  (expected: " << std::setw(10) << expected << ")\n";
    }

    // Orthogonality verification (numerical integration)
    std::cout << "\n--- Orthogonality Verification ---\n";
    std::cout << "(Integral should be 0 for m ≠ n)\n";
    int nPoints = 1000;
    Real dx = 2.0 / nPoints;
    
    for (int m = 0; m <= 2; ++m)
    {
        for (int n = m; n <= 2; ++n)
        {
            Real integral = 0.0;
            for (int i = 0; i < nPoints; ++i)
            {
                Real xi = -1.0 + (i + 0.5) * dx;
                integral += legendre.Evaluate(m, xi) * legendre.Evaluate(n, xi) * dx;
            }
            std::cout << "∫P_" << m << "(x)P_" << n << "(x)dx = " 
                      << std::setw(12) << std::setprecision(6) << integral;
            if (m == n)
                std::cout << "  (should be " << legendre.Normalization(n) << ")\n";
            else
                std::cout << "  (should be ~0)\n";
        }
    }

    std::cout << "\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
// Demo: Hermite Polynomials
// Domain: (-∞, ∞), Weight: w(x) = e^(-x²)
///////////////////////////////////////////////////////////////////////////////////////////
void Docs_Demo_HermiteBasis()
{
    std::cout << "========================================\n";
    std::cout << "   HERMITE POLYNOMIALS DEMO\n";
    std::cout << "========================================\n\n";

    HermiteBasis hermite;

    // Basic evaluation - first few Hermite polynomials  
    std::cout << "--- Polynomial Values at x = 1.0 ---\n";
    Real x = 1.0;
    for (int n = 0; n <= 5; ++n)
    {
        std::cout << "H_" << n << "(1.0) = " << std::setw(12) << std::setprecision(6)
                  << hermite.Evaluate(n, x) << "\n";
    }

    // Known values check
    std::cout << "\n--- Known Value Verification ---\n";
    std::cout << "H_0(x) = 1:       H_0(2) = " << hermite.Evaluate(0, 2.0) << "\n";
    std::cout << "H_1(x) = 2x:      H_1(2) = " << hermite.Evaluate(1, 2.0);
    std::cout << " (expected: " << 2.0*2.0 << ")\n";
    std::cout << "H_2(x) = 4x²-2:   H_2(2) = " << hermite.Evaluate(2, 2.0);
    std::cout << " (expected: " << 4.0*4.0 - 2.0 << ")\n";
    std::cout << "H_3(x) = 8x³-12x: H_3(2) = " << hermite.Evaluate(3, 2.0);
    std::cout << " (expected: " << 8.0*8.0 - 12.0*2.0 << ")\n";

    // Domain and weight
    std::cout << "\n--- Domain and Weight ---\n";
    std::cout << "Domain: (" << hermite.DomainMin() << ", " << hermite.DomainMax() << ")\n";
    std::cout << "Weight w(0) = " << hermite.WeightFunction(0.0) << " (should be 1.0)\n";
    std::cout << "Weight w(1) = " << hermite.WeightFunction(1.0) << " (should be e^-1 ≈ 0.368)\n";
    std::cout << "Weight w(2) = " << hermite.WeightFunction(2.0) << " (should be e^-4 ≈ 0.018)\n";

    // Normalization: ||Hₙ||² = 2ⁿn!√π
    std::cout << "\n--- Normalization ||Hₙ||² = 2ⁿn!√π ---\n";
    for (int n = 0; n <= 4; ++n)
    {
        // Calculate expected: 2^n * n! * sqrt(pi)
        Real factorial = 1.0;
        for (int k = 2; k <= n; ++k) factorial *= k;
        Real expected = std::pow(2.0, n) * factorial * std::sqrt(Constants::PI);
        
        std::cout << "||H_" << n << "||² = " << std::setw(12) << std::setprecision(6) 
                  << hermite.Normalization(n)
                  << "  (expected: " << std::setw(12) << expected << ")\n";
    }

    // Quantum harmonic oscillator wavefunctions
    std::cout << "\n--- Quantum Harmonic Oscillator Wavefunctions ---\n";
    std::cout << "ψₙ(x) ∝ Hₙ(x)e^(-x²/2)\n\n";
    std::cout << std::setw(6) << "x";
    for (int n = 0; n <= 3; ++n)
        std::cout << std::setw(12) << "ψ_" + std::to_string(n);
    std::cout << "\n";
    
    for (Real xi = -2.0; xi <= 2.01; xi += 0.5)
    {
        std::cout << std::setw(6) << std::setprecision(1) << std::fixed << xi;
        for (int n = 0; n <= 3; ++n)
        {
            Real psi = hermite.Evaluate(n, xi) * std::exp(-xi*xi/2.0);
            std::cout << std::setw(12) << std::setprecision(4) << psi;
        }
        std::cout << "\n";
    }

    std::cout << "\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
// Demo: Laguerre Polynomials
// Domain: [0, ∞), Weight: w(x) = e^(-x)
///////////////////////////////////////////////////////////////////////////////////////////
void Docs_Demo_LaguerreBasis()
{
    std::cout << "========================================\n";
    std::cout << "   LAGUERRE POLYNOMIALS DEMO\n";
    std::cout << "========================================\n\n";

    LaguerreBasis laguerre;

    // Basic evaluation - first few Laguerre polynomials
    std::cout << "--- Polynomial Values at x = 1.0 ---\n";
    Real x = 1.0;
    for (int n = 0; n <= 5; ++n)
    {
        std::cout << "L_" << n << "(1.0) = " << std::setw(12) << std::setprecision(6)
                  << laguerre.Evaluate(n, x) << "\n";
    }

    // Known values check
    std::cout << "\n--- Known Value Verification ---\n";
    std::cout << "L_0(x) = 1:       L_0(2) = " << laguerre.Evaluate(0, 2.0) << "\n";
    std::cout << "L_1(x) = 1-x:     L_1(2) = " << laguerre.Evaluate(1, 2.0);
    std::cout << " (expected: " << 1.0 - 2.0 << ")\n";
    std::cout << "L_2(x) = (x²-4x+2)/2: L_2(2) = " << laguerre.Evaluate(2, 2.0);
    std::cout << " (expected: " << (4.0 - 8.0 + 2.0)/2.0 << ")\n";
    
    // Property: Lₙ(0) = 1 for all n
    std::cout << "\n--- Property: Lₙ(0) = 1 ---\n";
    for (int n = 0; n <= 5; ++n)
    {
        std::cout << "L_" << n << "(0) = " << laguerre.Evaluate(n, 0.0) << "\n";
    }

    // Domain and weight
    std::cout << "\n--- Domain and Weight ---\n";
    std::cout << "Domain: [" << laguerre.DomainMin() << ", " << laguerre.DomainMax() << "]\n";
    std::cout << "Weight w(0) = " << laguerre.WeightFunction(0.0) << " (should be 1.0)\n";
    std::cout << "Weight w(1) = " << laguerre.WeightFunction(1.0) << " (should be e^-1 ≈ 0.368)\n";
    std::cout << "Weight w(2) = " << laguerre.WeightFunction(2.0) << " (should be e^-2 ≈ 0.135)\n";

    // Normalization: ||Lₙ||² = 1
    std::cout << "\n--- Normalization ||Lₙ||² = 1 ---\n";
    for (int n = 0; n <= 4; ++n)
    {
        std::cout << "||L_" << n << "||² = " << laguerre.Normalization(n) << "\n";
    }

    // Hydrogen atom radial wavefunctions
    std::cout << "\n--- Hydrogen Atom Radial Wavefunctions ---\n";
    std::cout << "Radial part involves Lₙ(2r/na₀)·e^(-r/na₀)\n";
    std::cout << "Showing e^(-x)·Lₙ(x) for first few n:\n\n";
    
    std::cout << std::setw(6) << "x";
    for (int n = 0; n <= 3; ++n)
        std::cout << std::setw(12) << "e^-x·L_" + std::to_string(n);
    std::cout << "\n";
    
    for (Real xi = 0.0; xi <= 5.01; xi += 0.5)
    {
        std::cout << std::setw(6) << std::setprecision(1) << std::fixed << xi;
        for (int n = 0; n <= 3; ++n)
        {
            Real val = std::exp(-xi) * laguerre.Evaluate(n, xi);
            std::cout << std::setw(12) << std::setprecision(4) << val;
        }
        std::cout << "\n";
    }

    std::cout << "\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
// Demo: Chebyshev Polynomials (First Kind)
// Domain: [-1, 1], Weight: w(x) = 1/√(1-x²)
///////////////////////////////////////////////////////////////////////////////////////////
void Docs_Demo_ChebyshevBasis()
{
    std::cout << "========================================\n";
    std::cout << "   CHEBYSHEV POLYNOMIALS DEMO\n";
    std::cout << "========================================\n\n";

    ChebyshevBasis chebyshev;

    // Basic evaluation - first few Chebyshev polynomials
    std::cout << "--- Polynomial Values at x = 0.5 ---\n";
    Real x = 0.5;
    for (int n = 0; n <= 5; ++n)
    {
        std::cout << "T_" << n << "(0.5) = " << std::setw(12) << std::setprecision(8)
                  << chebyshev.Evaluate(n, x) << "\n";
    }

    // Special property: Tₙ(cos θ) = cos(nθ)
    std::cout << "\n--- Special Property: Tₙ(cos θ) = cos(nθ) ---\n";
    Real theta = Constants::PI / 4.0;  // 45 degrees
    Real cos_theta = std::cos(theta);
    std::cout << "θ = π/4, cos(θ) = " << cos_theta << "\n\n";
    
    for (int n = 0; n <= 4; ++n)
    {
        Real Tn = chebyshev.Evaluate(n, cos_theta);
        Real cos_n_theta = std::cos(n * theta);
        std::cout << "T_" << n << "(cos(π/4)) = " << std::setw(10) << std::setprecision(6) << Tn
                  << "   cos(" << n << "·π/4) = " << std::setw(10) << cos_n_theta << "\n";
    }

    // Bounded property: |Tₙ(x)| ≤ 1 for x ∈ [-1, 1]
    std::cout << "\n--- Bounded Property: |Tₙ(x)| ≤ 1 ---\n";
    std::cout << "Evaluating T_5(x) across [-1, 1]:\n";
    std::cout << std::setw(8) << "x" << std::setw(12) << "T_5(x)" << "\n";
    for (Real xi = -1.0; xi <= 1.01; xi += 0.25)
    {
        std::cout << std::setw(8) << std::setprecision(2) << std::fixed << xi
                  << std::setw(12) << std::setprecision(6) << chebyshev.Evaluate(5, xi) << "\n";
    }

    // Domain and weight
    std::cout << "\n--- Domain and Weight ---\n";
    std::cout << "Domain: [" << chebyshev.DomainMin() << ", " << chebyshev.DomainMax() << "]\n";
    std::cout << "Weight w(0) = " << chebyshev.WeightFunction(0.0) << " (1/√1 = 1.0)\n";
    std::cout << "Weight w(0.5) = " << chebyshev.WeightFunction(0.5) 
              << " (1/√0.75 ≈ 1.155)\n";
    std::cout << "Weight w(0.9) = " << chebyshev.WeightFunction(0.9) 
              << " (1/√0.19 ≈ 2.294)\n";

    // Normalization: ||T₀||² = π, ||Tₙ||² = π/2 for n > 0
    std::cout << "\n--- Normalization ---\n";
    std::cout << "||T_0||² = " << chebyshev.Normalization(0) << " (should be π ≈ 3.14159)\n";
    for (int n = 1; n <= 4; ++n)
    {
        std::cout << "||T_" << n << "||² = " << chebyshev.Normalization(n) 
                  << " (should be π/2 ≈ 1.5708)\n";
    }

    // Zeros and extrema
    std::cout << "\n--- Zeros and Extrema of T_4(x) ---\n";
    std::cout << "Zeros: xₖ = cos((2k-1)π/(2n)), k = 1,...,n\n";
    int n = 4;
    std::cout << "Zeros of T_4:\n";
    for (int k = 1; k <= n; ++k)
    {
        Real zero = std::cos((2.0 * k - 1.0) * Constants::PI / (2.0 * n));
        std::cout << "  x_" << k << " = " << std::setw(10) << std::setprecision(6) << zero
                  << "  T_4(" << std::setprecision(4) << zero << ") = " 
                  << std::setprecision(2) << std::scientific << chebyshev.Evaluate(4, zero) << "\n";
    }
    
    std::cout << "\nExtrema: xₖ = cos(kπ/n), k = 0,...,n\n";
    std::cout << "Extrema of T_4 (values ±1):\n";
    for (int k = 0; k <= n; ++k)
    {
        Real extrema = std::cos(static_cast<Real>(k) * Constants::PI / n);
        std::cout << "  x_" << k << " = " << std::setw(10) << std::setprecision(6) << std::fixed << extrema
                  << "  T_4(" << std::setprecision(4) << extrema << ") = " 
                  << std::setw(6) << std::setprecision(2) << chebyshev.Evaluate(4, extrema) << "\n";
    }

    std::cout << "\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
// Demo: Function Expansion Using Orthogonal Bases
///////////////////////////////////////////////////////////////////////////////////////////
void Docs_Demo_FunctionExpansion()
{
    std::cout << "========================================\n";
    std::cout << "   FUNCTION EXPANSION DEMO\n";
    std::cout << "========================================\n\n";

    // Expand f(x) = x² using Legendre polynomials
    std::cout << "--- Expanding f(x) = x² in Legendre Polynomials ---\n";
    std::cout << "f(x) = Σ cₙPₙ(x) where cₙ = <f, Pₙ>/||Pₙ||²\n\n";

    LegendreBasis legendre;
    
    // f(x) = x²
    std::function<Real(Real)> f = [](Real x) { return x * x; };
    RealFunctionFromStdFunc fWrapper(f);

    // Compute coefficients using numerical integration
    int nPoints = 1000;
    Real dx = 2.0 / nPoints;
    
    std::vector<Real> coeffs(5);
    for (int n = 0; n < 5; ++n)
    {
        Real inner_product = 0.0;
        for (int i = 0; i < nPoints; ++i)
        {
            Real xi = -1.0 + (i + 0.5) * dx;
            inner_product += f(xi) * legendre.Evaluate(n, xi) * dx;
        }
        coeffs[n] = inner_product / legendre.Normalization(n);
    }

    std::cout << "Expansion coefficients:\n";
    for (int n = 0; n < 5; ++n)
    {
        std::cout << "c_" << n << " = " << std::setw(12) << std::setprecision(6) << coeffs[n] << "\n";
    }

    std::cout << "\nTheoretical: x² = (1/3)P₀(x) + 0·P₁(x) + (2/3)P₂(x)\n";
    std::cout << "c_0 = 1/3 ≈ 0.333333, c_1 = 0, c_2 = 2/3 ≈ 0.666667\n";

    // Verify reconstruction
    std::cout << "\n--- Reconstruction Verification ---\n";
    std::cout << std::setw(8) << "x" << std::setw(12) << "f(x)=x²" 
              << std::setw(14) << "Σcₙ·Pₙ(x)" << std::setw(12) << "Error\n";
    
    for (Real xi = -1.0; xi <= 1.01; xi += 0.25)
    {
        Real f_exact = xi * xi;
        Real f_approx = 0.0;
        for (int n = 0; n < 5; ++n)
        {
            f_approx += coeffs[n] * legendre.Evaluate(n, xi);
        }
        std::cout << std::setw(8) << std::setprecision(2) << std::fixed << xi
                  << std::setw(12) << std::setprecision(4) << f_exact
                  << std::setw(14) << std::setprecision(6) << f_approx
                  << std::setw(12) << std::scientific << std::setprecision(2) 
                  << std::abs(f_exact - f_approx) << "\n";
    }

    std::cout << "\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
// Demo: Orthogonality Comparison
///////////////////////////////////////////////////////////////////////////////////////////
void Docs_Demo_OrthogonalityComparison()
{
    std::cout << "========================================\n";
    std::cout << "   ORTHOGONALITY COMPARISON\n";
    std::cout << "========================================\n\n";

    std::cout << "Each basis is orthogonal under its specific weight function.\n\n";

    std::cout << "┌─────────────┬──────────────────┬───────────────────┬────────────────┐\n";
    std::cout << "│ Basis       │ Domain           │ Weight w(x)       │ Inner Product  │\n";
    std::cout << "├─────────────┼──────────────────┼───────────────────┼────────────────┤\n";
    std::cout << "│ Legendre    │ [-1, 1]          │ 1                 │ ∫₋₁¹ f·g dx    │\n";
    std::cout << "│ Hermite     │ (-∞, ∞)          │ e^(-x²)           │ ∫ f·g·w dx     │\n";
    std::cout << "│ Laguerre    │ [0, ∞)           │ e^(-x)            │ ∫₀^∞ f·g·w dx  │\n";
    std::cout << "│ Chebyshev   │ [-1, 1]          │ 1/√(1-x²)         │ ∫₋₁¹ f·g·w dx  │\n";
    std::cout << "└─────────────┴──────────────────┴───────────────────┴────────────────┘\n";

    std::cout << "\n--- Physical Applications ---\n";
    std::cout << "• Legendre:  Potential theory, multipole expansion, angular momentum\n";
    std::cout << "• Hermite:   Quantum harmonic oscillator, probability distributions\n";
    std::cout << "• Laguerre:  Hydrogen atom wavefunctions, heat conduction\n";
    std::cout << "• Chebyshev: Polynomial approximation (minimax), spectral methods\n";

    std::cout << "\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
// Main Demo Entry Point
///////////////////////////////////////////////////////////////////////////////////////////
void Docs_Demo_OrthogonalBases()
{
    std::cout << "\n";
    std::cout << "╔════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║         ORTHOGONAL POLYNOMIAL BASES - DOCUMENTATION DEMOS      ║\n";
    std::cout << "║                                                                ║\n";
    std::cout << "║  Demonstrates:                                                 ║\n";
    std::cout << "║    • Legendre polynomials Pₙ(x)                                ║\n";
    std::cout << "║    • Hermite polynomials Hₙ(x)                                 ║\n";
    std::cout << "║    • Laguerre polynomials Lₙ(x)                                ║\n";
    std::cout << "║    • Chebyshev polynomials Tₙ(x)                               ║\n";
    std::cout << "║    • Function expansion using orthogonal bases                 ║\n";
    std::cout << "╚════════════════════════════════════════════════════════════════╝\n";
    std::cout << "\n";

    Docs_Demo_LegendreBasis();
    Docs_Demo_HermiteBasis();
    Docs_Demo_LaguerreBasis();
    Docs_Demo_ChebyshevBasis();
    Docs_Demo_FunctionExpansion();
    Docs_Demo_OrthogonalityComparison();

    std::cout << "========================================\n";
    std::cout << "   ORTHOGONAL BASES DEMOS COMPLETE\n";
    std::cout << "========================================\n\n";
}
