///////////////////////////////////////////////////////////////////////////////////////////
// README Example 0: Gauss's Divergence Theorem Verification
// 
// This example demonstrates MML's capabilities by verifying one of the fundamental
// theorems of vector calculus: Gauss's Divergence Theorem
//
//    ∫∫∫(∇·F)dV = ∮∮(F·n̂)dS
//
// The volume integral of divergence equals the surface integral of flux through
// a closed surface.
///////////////////////////////////////////////////////////////////////////////////////////

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/Vector.h"
#include "base/VectorN.h"
#include "base/Function.h"
#include "base/Geometry3D.h"
#include "base/Geometry3DBodies.h"

#include "core/Integration.h"
#include "core/Integration/SurfaceIntegration.h"
#include "core/FieldOperations.h"
#endif

#include <iomanip>

using namespace MML;

void Readme_demonstration()
{
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****     README - Gauss's Divergence Theorem Demonstration        ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;
    std::cout << std::endl;

    // ============================================================================
    // Verify Gauss's Divergence Theorem: ∫∫∫(∇·F)dV = ∮∮(F·n̂)dS
    // ============================================================================

    // Define a 3D vector field F(x,y,z) = (x², xy, z)
    VectorFunction<3> F([](const VectorN<Real, 3>& p) {
        return VectorN<Real, 3>{ p[0]*p[0], p[0]*p[1], p[2] };
    });

    // Compute divergence analytically:
    // ∇·F = ∂(x²)/∂x + ∂(xy)/∂y + ∂z/∂z = 2x + x + 1 = 3x + 1
    ScalarFunction<3> divF([](const VectorN<Real, 3>& p) { 
        return 3*p[0] + 1; 
    });

    // Define limits for unit cube [0,1]³
    auto y_low = [](Real x) { return 0.0; };
    auto y_high = [](Real x) { return 1.0; };
    auto z_low = [](Real x, Real y) { return 0.0; };
    auto z_high = [](Real x, Real y) { return 1.0; };

    // Integrate divergence over unit cube [0,1]³
    // ∫∫∫(3x + 1)dV = ∫₀¹∫₀¹∫₀¹(3x + 1)dzdydx
    //               = ∫₀¹∫₀¹(3x + 1)dydx
    //               = ∫₀¹(3x + 1)dx
    //               = [3x²/2 + x]₀¹ = 3/2 + 1 = 2.5
    Real volIntegral = Integrate3D(divF, 0.0, 1.0, y_low, y_high, z_low, z_high);

    // Surface integral over 6 faces of the unit cube (flux through closed surface)
    // Unit cube centered at (0.5, 0.5, 0.5) with side length 1
    Cube3D unitCube(1.0, Point3Cartesian(0.5, 0.5, 0.5));
    Real fluxIntegral = SurfaceIntegration::SurfaceIntegral(F, unitCube, 1e-6);

    std::cout << "Vector field F(x,y,z) = (x², xy, z)" << std::endl;
    std::cout << "Divergence ∇·F = 2x + x + 1 = 3x + 1" << std::endl;
    std::cout << "Integration region: unit cube [0,1]³" << std::endl;
    std::cout << std::endl;

    std::cout << std::setprecision(15);
    std::cout << "Volume integral   ∫∫∫(∇·F)dV = " << volIntegral << std::endl;
    std::cout << "Surface integral  ∮∮(F·n̂)dS  = " << fluxIntegral << std::endl;
    std::cout << "(Gaussian quadrature achieves exact results for low-degree polynomials)" << std::endl;
    std::cout << std::endl;

    // ============================================================================
    // Second example with NON-separable trigonometric field
    // ============================================================================
    std::cout << "--- Example 2: Non-separable Trigonometric Field ---" << std::endl;
    
    // F(x,y,z) = (sin(xy), cos(yz), exp(xz))
    // ∇·F = y*cos(xy) - z*sin(yz) + x*exp(xz)
    VectorFunction<3> F2([](const VectorN<Real, 3>& p) {
        return VectorN<Real, 3>{ sin(p[0]*p[1]), cos(p[1]*p[2]), exp(p[0]*p[2]) };
    });
    
    ScalarFunction<3> divF2([](const VectorN<Real, 3>& p) { 
        return p[1]*cos(p[0]*p[1]) - p[2]*sin(p[1]*p[2]) + p[0]*exp(p[0]*p[2]); 
    });
    
    IntegrationResult volResult2 = Integrate3D(divF2, 0.0, 1.0, y_low, y_high, z_low, z_high);
    Real fluxIntegral2 = SurfaceIntegration::SurfaceIntegral(F2, unitCube, 1e-8);
    
    std::cout << std::setprecision(15);
    std::cout << "Vector field F(x,y,z) = (sin(xy), cos(yz), exp(xz))" << std::endl;
    std::cout << "Volume integral   = " << volResult2.value << std::endl;
    std::cout << "Surface integral  = " << fluxIntegral2 << std::endl;
    std::cout << "Difference        = " << std::abs(volResult2.value - fluxIntegral2) << std::endl;
    std::cout << std::endl;

    // ============================================================================
    // Third example: NUMERICAL divergence (no analytical formula needed!)
    // ============================================================================
    std::cout << "--- Example 3: Numerically Computed Divergence ---" << std::endl;
    std::cout << "Same field F(x,y,z) = (sin(xy), cos(yz), exp(xz))" << std::endl;
    std::cout << "But divergence computed NUMERICALLY from the vector field!" << std::endl;
    std::cout << std::endl;
    
    // Count function evaluations to show the computational cost
    int analyticalDivEvals = 0;
    int numericalDivEvals = 0;
    int vectorFieldEvals = 0;
    
    // Analytical divergence with evaluation counter
    ScalarFunctionFromStdFunc<3> divF2_counted([&analyticalDivEvals](const VectorN<Real, 3>& p) {
        analyticalDivEvals++;
        return p[1]*cos(p[0]*p[1]) - p[2]*sin(p[1]*p[2]) + p[0]*exp(p[0]*p[2]); 
    });
    
    // Numerical divergence with evaluation counters
    // Note: DivCart computes 6 partial derivatives, each requiring 2 field evaluations
    VectorFunctionFromStdFunc<3> F2_counted([&vectorFieldEvals](const VectorN<Real, 3>& p) {
        vectorFieldEvals++;
        return VectorN<Real, 3>{ sin(p[0]*p[1]), cos(p[1]*p[2]), exp(p[0]*p[2]) };
    });
    
    ScalarFunctionFromStdFunc<3> divF2_numerical([&F2_counted, &numericalDivEvals](const VectorN<Real, 3>& p) {
        numericalDivEvals++;
        return VectorFieldOperations::DivCart<3>(F2_counted, p);
    });
    
    IntegrationResult volResultAnalytical = Integrate3D(divF2_counted, 0.0, 1.0, y_low, y_high, z_low, z_high);
    IntegrationResult volResultNumerical = Integrate3D(divF2_numerical, 0.0, 1.0, y_low, y_high, z_low, z_high);
    
    std::cout << "Volume integral (numerical div)  = " << volResultNumerical.value << std::endl;
    std::cout << "Volume integral (analytical div) = " << volResultAnalytical.value << std::endl;
    std::cout << "Surface integral                 = " << fluxIntegral2 << std::endl;
    std::cout << std::endl;
    
    // Show function evaluation counts
    std::cout << "Function evaluations:" << std::endl;
    std::cout << "  Analytical divergence evals:   " << analyticalDivEvals << std::endl;
    std::cout << "  Numerical divergence evals:    " << numericalDivEvals << std::endl;
    std::cout << "  Vector field evals (for div):  " << vectorFieldEvals << std::endl;
    std::cout << "  (Gauss-10 uses 10³ = 1000 quadrature points)" << std::endl;
    std::cout << "  (Numerical divergence needs ~6× more field evals per point)" << std::endl;
    std::cout << std::endl;
    
    std::cout << "Difference (numerical vs analytical div): " << std::abs(volResultNumerical.value - volResultAnalytical.value) << std::endl;
    std::cout << "Difference (numerical div vs surface):    " << std::abs(volResultNumerical.value - fluxIntegral2) << std::endl;
    std::cout << std::endl;

    // ============================================================================
    // Fourth example: Compare different integration methods
    // ============================================================================
    std::cout << "--- Example 4: Comparing Integration Methods ---" << std::endl;
    std::cout << "Same numerical divergence, different integration methods:" << std::endl;
    std::cout << std::endl;
    
    // Reset counter for fresh measurement
    VectorFunctionFromStdFunc<3> F2_for_methods([](const VectorN<Real, 3>& p) {
        return VectorN<Real, 3>{ sin(p[0]*p[1]), cos(p[1]*p[2]), exp(p[0]*p[2]) };
    });
    ScalarFunctionFromStdFunc<3> divF2_for_methods([&F2_for_methods](const VectorN<Real, 3>& p) {
        return VectorFieldOperations::DivCart<3>(F2_for_methods, p);
    });
    
    IntegrationResult resultGauss = Integrate3D(divF2_for_methods, GAUSS10, 0.0, 1.0, y_low, y_high, z_low, z_high);
    IntegrationResult resultSimpson = Integrate3D(divF2_for_methods, SIMPSON, 0.0, 1.0, y_low, y_high, z_low, z_high);
    IntegrationResult resultTrap = Integrate3D(divF2_for_methods, TRAP, 0.0, 1.0, y_low, y_high, z_low, z_high);
    IntegrationResult resultRomberg = Integrate3D(divF2_for_methods, ROMBERG, 0.0, 1.0, y_low, y_high, z_low, z_high);
    
    std::cout << std::setprecision(14);
    std::cout << "  GAUSS10: " << resultGauss.value << "  (error vs surface: " << std::abs(resultGauss.value - fluxIntegral2) << ")" << std::endl;
    std::cout << "  SIMPSON: " << resultSimpson.value << "  (error vs surface: " << std::abs(resultSimpson.value - fluxIntegral2) << ")" << std::endl;
    std::cout << "  ROMBERG: " << resultRomberg.value << "  (error vs surface: " << std::abs(resultRomberg.value - fluxIntegral2) << ")" << std::endl;
    std::cout << "  TRAP:    " << resultTrap.value << "  (error vs surface: " << std::abs(resultTrap.value - fluxIntegral2) << ")" << std::endl;
    std::cout << std::endl;

    Real error = std::abs(volIntegral - fluxIntegral);
    bool verified = error < 1e-4;
    
    std::cout << "Analytical result: 2.5" << std::endl;
    std::cout << "Difference: |vol - surf| = " << error << std::endl;
    std::cout << "Gauss's Theorem verified: " << (verified ? "YES ✓" : "NO ✗") << std::endl;
    std::cout << std::endl;

    // ============================================================================
    // Additional demonstration: Verify with analytical derivation
    // ============================================================================
    std::cout << "--- Analytical Verification ---" << std::endl;
    std::cout << "For F = (x², xy, z) over [0,1]³:" << std::endl;
    std::cout << std::endl;
    std::cout << "Volume integral of divergence:" << std::endl;
    std::cout << "  ∫₀¹∫₀¹∫₀¹(3x+1) dz dy dx = ∫₀¹(3x+1) dx = 3/2 + 1 = 2.5" << std::endl;
    std::cout << std::endl;
    std::cout << "Surface integral (6 faces):" << std::endl;
    std::cout << "  Face x=0: ∫∫(0,0,z)·(-1,0,0) dydz = 0" << std::endl;
    std::cout << "  Face x=1: ∫∫(1,y,z)·(1,0,0) dydz = ∫₀¹∫₀¹ 1 dydz = 1" << std::endl;
    std::cout << "  Face y=0: ∫∫(x²,0,z)·(0,-1,0) dxdz = 0" << std::endl;
    std::cout << "  Face y=1: ∫∫(x²,x,z)·(0,1,0) dxdz = ∫₀¹∫₀¹ x dxdz = 0.5" << std::endl;
    std::cout << "  Face z=0: ∫∫(x²,xy,0)·(0,0,-1) dxdy = 0" << std::endl;
    std::cout << "  Face z=1: ∫∫(x²,xy,1)·(0,0,1) dxdy = ∫₀¹∫₀¹ 1 dxdy = 1" << std::endl;
    std::cout << "  Total flux = 0 + 1 + 0 + 0.5 + 0 + 1 = 2.5 ✓" << std::endl;
    std::cout << std::endl;
}
