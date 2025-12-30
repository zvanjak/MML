#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/DiracDeltaFunction.h"
#include "base/Function.h"
#include "core/Integration.h"
#endif

#include <iostream>
#include <iomanip>

using namespace MML;

///////////////////////////////////////////////////////////////////////////////
/// Dirac Delta Function Approximations Demo
/// 
/// Demonstrates the four Dirac delta approximations:
/// - DiracStep: Rectangular (box) approximation
/// - DiracExp: Gaussian approximation  
/// - DiracSqr: Lorentzian (Cauchy) approximation
/// - DiracSin: Sinc approximation
///////////////////////////////////////////////////////////////////////////////

void Docs_Demo_DiracStep()
{
    std::cout << "\n=== DiracStep - Rectangular Approximation ===" << std::endl;
    std::cout << "Formula: delta(x) = N for |x| < 1/(2N), else 0" << std::endl;
    std::cout << "Properties: Width = 1/N, Height = N, Area = 1 (exact)" << std::endl;
    
    DiracStep delta10(10);
    DiracStep delta100(100);
    DiracStep delta1000(1000);
    
    std::cout << "\nN = 10:   delta(0) = " << delta10(0.0) 
              << ", delta(0.1) = " << delta10(0.1)
              << ", width = " << 1.0/10 << std::endl;
    
    std::cout << "N = 100:  delta(0) = " << delta100(0.0) 
              << ", delta(0.01) = " << delta100(0.01)
              << ", width = " << 1.0/100 << std::endl;
              
    std::cout << "N = 1000: delta(0) = " << delta1000(0.0) 
              << ", delta(0.001) = " << delta1000(0.001)
              << ", width = " << 1.0/1000 << std::endl;
    
    // Test boundary behavior
    DiracStep delta(100);
    Real halfWidth = 1.0 / (2 * 100);  // 0.005
    std::cout << "\nBoundary test (N=100, half-width=" << halfWidth << "):" << std::endl;
    std::cout << "  delta(" << halfWidth - 0.001 << ") = " << delta(halfWidth - 0.001) << " (inside)" << std::endl;
    std::cout << "  delta(" << halfWidth + 0.001 << ") = " << delta(halfWidth + 0.001) << " (outside)" << std::endl;
}

void Docs_Demo_DiracExp()
{
    std::cout << "\n=== DiracExp - Gaussian Approximation ===" << std::endl;
    std::cout << "Formula: delta(x) = N/sqrt(2*pi) * exp(-x^2 * N^2)" << std::endl;
    std::cout << "Properties: Peak = N/sqrt(2*pi), Std Dev = 1/N" << std::endl;
    
    DiracExp delta10(10);
    DiracExp delta100(100);
    DiracExp delta1000(1000);
    
    Real expectedPeak10 = 10.0 / std::sqrt(2 * Constants::PI);
    Real expectedPeak100 = 100.0 / std::sqrt(2 * Constants::PI);
    Real expectedPeak1000 = 1000.0 / std::sqrt(2 * Constants::PI);
    
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "\nN = 10:   delta(0) = " << delta10(0.0) 
              << " (expected: " << expectedPeak10 << ")" << std::endl;
    std::cout << "N = 100:  delta(0) = " << delta100(0.0) 
              << " (expected: " << expectedPeak100 << ")" << std::endl;
    std::cout << "N = 1000: delta(0) = " << delta1000(0.0) 
              << " (expected: " << expectedPeak1000 << ")" << std::endl;
    
    // Show decay
    DiracExp delta(100);
    std::cout << "\nGaussian decay (N=100):" << std::endl;
    for (Real x = 0.0; x <= 0.05; x += 0.01) {
        std::cout << "  delta(" << x << ") = " << delta(x) << std::endl;
    }
}

void Docs_Demo_DiracSqr()
{
    std::cout << "\n=== DiracSqr - Lorentzian (Cauchy) Approximation ===" << std::endl;
    std::cout << "Formula: delta(x) = N / (pi * (1 + N^2 * x^2))" << std::endl;
    std::cout << "Properties: Peak = N/pi, HWHM = 1/N" << std::endl;
    
    DiracSqr delta10(10);
    DiracSqr delta100(100);
    DiracSqr delta1000(1000);
    
    Real expectedPeak10 = 10.0 / Constants::PI;
    Real expectedPeak100 = 100.0 / Constants::PI;
    Real expectedPeak1000 = 1000.0 / Constants::PI;
    
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "\nN = 10:   delta(0) = " << delta10(0.0) 
              << " (expected: " << expectedPeak10 << ")" << std::endl;
    std::cout << "N = 100:  delta(0) = " << delta100(0.0) 
              << " (expected: " << expectedPeak100 << ")" << std::endl;
    std::cout << "N = 1000: delta(0) = " << delta1000(0.0) 
              << " (expected: " << expectedPeak1000 << ")" << std::endl;
    
    // Show heavier tails compared to Gaussian
    DiracSqr deltaL(100);
    DiracExp deltaG(100);
    std::cout << "\nLorentzian vs Gaussian tails (N=100):" << std::endl;
    std::cout << "  x\t\tLorentzian\tGaussian" << std::endl;
    for (Real x = 0.0; x <= 0.1; x += 0.02) {
        std::cout << "  " << x << "\t\t" << deltaL(x) << "\t\t" << deltaG(x) << std::endl;
    }
}

void Docs_Demo_DiracSin()
{
    std::cout << "\n=== DiracSin - Sinc Approximation ===" << std::endl;
    std::cout << "Formula: delta(x) = sin(N*x) / (pi * x)" << std::endl;
    std::cout << "Properties: Oscillating tails, limit at x=0 is N/pi" << std::endl;
    
    DiracSin delta10(10);
    DiracSin delta100(100);
    
    // Note: undefined at x=0, show values near zero
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "\nN = 10 (limit at 0 = " << 10.0/Constants::PI << "):" << std::endl;
    for (Real x = 0.001; x <= 0.5; x += 0.1) {
        std::cout << "  delta(" << x << ") = " << delta10(x) << std::endl;
    }
    
    // Show oscillations
    std::cout << "\nOscillating behavior (N=10):" << std::endl;
    for (Real x = 0.1; x <= 1.0; x += 0.1) {
        std::cout << "  delta(" << x << ") = " << delta10(x) << std::endl;
    }
}

void Docs_Demo_DiracDelta_Integration()
{
    std::cout << "\n=== Integration with Dirac Delta ===" << std::endl;
    std::cout << "Property: integral[f(x) * delta(x-a)] = f(a)" << std::endl;
    
    // Test function: f(x) = x^2
    auto testFunc = [](Real x) { return x * x; };
    
    // Sample f at x = 2 using delta approximation
    Real a = 2.0;
    Real expected = testFunc(a);  // f(2) = 4
    
    std::cout << "\nSampling f(x) = x^2 at x = " << a << " (expected: " << expected << ")" << std::endl;
    
    // Try different N values with DiracExp
    for (int N = 10; N <= 1000; N *= 10) {
        DiracExp delta(N);
        
        // Create integrand: f(x) * delta(x - a)
        RealFunctionFromStdFunc integrand([&testFunc, &delta, a](Real x) {
            return testFunc(x) * delta(x - a);
        });
        
        // Integrate over a wide range centered at a
        Real result = IntegrateTrap(integrand, a - 5.0, a + 5.0, 10000);
        Real error = std::abs(result - expected);
        
        std::cout << "  N = " << std::setw(4) << N 
                  << ": integral = " << std::setw(10) << result 
                  << ", error = " << error << std::endl;
    }
}

void Docs_Demo_DiracDelta_Comparison()
{
    std::cout << "\n=== Comparison of All Approximations ===" << std::endl;
    std::cout << "All with N = 100" << std::endl;
    
    DiracStep deltaStep(100);
    DiracExp deltaExp(100);
    DiracSqr deltaSqr(100);
    DiracSin deltaSin(100);
    
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "\n  x\t\tStep\t\tExp\t\tSqr\t\tSin" << std::endl;
    std::cout << "  ----\t\t----\t\t---\t\t---\t\t---" << std::endl;
    
    std::vector<Real> xValues = {0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1};
    for (Real x : xValues) {
        std::cout << "  " << x << "\t\t" 
                  << deltaStep(x) << "\t" 
                  << deltaExp(x) << "\t" 
                  << deltaSqr(x) << "\t"
                  << deltaSin(x) << std::endl;
    }
    
    // Peak values at x = 0
    std::cout << "\nPeak values at x = 0:" << std::endl;
    std::cout << "  Step: " << deltaStep(0.0) << " (exact: 100)" << std::endl;
    std::cout << "  Exp:  " << deltaExp(0.0) << " (exact: " << 100.0/std::sqrt(2*Constants::PI) << ")" << std::endl;
    std::cout << "  Sqr:  " << deltaSqr(0.0) << " (exact: " << 100.0/Constants::PI << ")" << std::endl;
    std::cout << "  Sin:  undefined (limit: " << 100.0/Constants::PI << ")" << std::endl;
}

void Docs_Demo_DiracDelta()
{
    std::cout << "**********************************************************************" << std::endl;
    std::cout << "***    Docs Demo - Dirac Delta Function Approximations             ***" << std::endl;
    std::cout << "**********************************************************************" << std::endl;
    
    Docs_Demo_DiracStep();
    Docs_Demo_DiracExp();
    Docs_Demo_DiracSqr();
    Docs_Demo_DiracSin();
    Docs_Demo_DiracDelta_Comparison();
    Docs_Demo_DiracDelta_Integration();
    
    std::cout << "\n**********************************************************************" << std::endl;
}
