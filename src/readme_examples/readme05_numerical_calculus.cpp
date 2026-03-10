///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        readme05_numerical_calculus.cpp                                     ///
///  Description: README Numerical Calculus section demo                              ///
///               Demonstrates derivatives and integration                            ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///////////////////////////////////////////////////////////////////////////////////////////

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/Vector/Vector.h"
#include "base/Vector/VectorN.h"
#include "base/Function.h"
#include "core/Derivation.h"
#include "core/Integration.h"
#include "core/Integration/MonteCarloIntegration.h"
#endif

#include <iomanip>

using namespace MML;

void Readme_NumericalCalculus()
{
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****            README - Numerical Calculus                        ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    // Compare derivative orders on f(x) = sin(x) at x = 1.0
    // Analytical derivative: f'(1) = cos(1) ≈ 0.5403023058681398
    RealFunction f{[](Real x) { return std::sin(x); }};
    Real x = 1.0;
    Real analytical = std::cos(1.0);

    Real der1 = Derivation::NDer1(f, x);   // O(h) - forward difference
    Real der2 = Derivation::NDer2(f, x);   // O(h²) - central difference
    Real der4 = Derivation::NDer4(f, x);   // O(h⁴) - 5-point stencil
    Real der6 = Derivation::NDer6(f, x);   // O(h⁶) - 7-point stencil
    Real der8 = Derivation::NDer8(f, x);   // O(h⁸) - 9-point stencil

    std::cout << std::scientific << std::setprecision(10);
    std::cout << "NDer1 error: " << std::abs(der1 - analytical) << std::endl;
    std::cout << "NDer2 error: " << std::abs(der2 - analytical) << std::endl;
    std::cout << "NDer8 error: " << std::abs(der8 - analytical) << std::endl;

    // Numerical Integration
    // Integrate f(x) = sin(x) from 0 to π — Analytical result: 2.0
    Real a = 0.0, b = Constants::PI;

    auto trap = IntegrateTrap(f, a, b, 1e-8);
    auto simp = IntegrateSimpson(f, a, b, 1e-8);
    auto romb = IntegrateRomberg(f, a, b, 1e-10);

    std::cout << std::fixed << std::setprecision(10);
    std::cout << "Trapezoid: " << trap.value << " (error: " << trap.error_estimate << ")\n";
    std::cout << "Simpson:   " << simp.value << " (error: " << simp.error_estimate << ")\n";
    std::cout << "Romberg:   " << romb.value << " (error: " << romb.error_estimate << ")\n";

    // Monte Carlo for high-dimensional integrals
    ScalarFunction<3> volume_func([](const VectorN<Real, 3>& v) {
        return v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
    });

    MonteCarloIntegrator<3> mc_integrator;
    VectorN<Real, 3> lower{0, 0, 0}, upper{1, 1, 1};
    auto mc_result = mc_integrator.integrate(volume_func, lower, upper, 
                                              MonteCarloConfig().samples(100000));
    std::cout << "MC estimate: " << mc_result.value 
              << " +/- " << mc_result.error_estimate << std::endl;

/* Expected OUTPUT (approx):
    NDer1 error: 1.4355650446e-08
    NDer2 error: 6.2718719107e-12
    NDer8 error: 3.9968028887e-15
    Trapezoid: 1.9999999939 (error: 0.0000000184)
    Simpson:   2.0000000003 (error: 0.0000000038)
    Romberg:   2.0000000000 (error: 0.0000000001)
    MC estimate: ~1.0 +/- ~0.001 (analytical: 1.0)
*/
}
