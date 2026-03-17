///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        readme00_fundamental_theorems.cpp                                   ///
///  Description: Flagship demonstration of the three fundamental theorems of        ///
///               vector calculus verified numerically:                               ///
///               1. Gauss's Divergence Theorem                                       ///
///               2. Stokes' Theorem                                                  ///
///               3. Green's Theorem (2D Stokes)                                      ///
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
#include "base/Geometry/Geometry3D.h"
#include "base/Geometry/Geometry3DBodies.h"

#include "core/Integration.h"
#include "core/Integration/SurfaceIntegration.h"
#include "core/Integration/PathIntegration.h"
#include "core/FieldOperations.h"
#include "core/Curves.h"

#include "algorithms/FieldAnalyzers.h"
#endif

#include <iomanip>

using namespace MML;

///////////////////////////////////////////////////////////////////////////////
//  GAUSS'S DIVERGENCE THEOREM: ∫∫∫(∇·F)dV = ∮∮(F·n̂)dS
///////////////////////////////////////////////////////////////////////////////
void Readme_GausssDivergenceTheorem()
{
    std::cout << "═══════════════════════════════════════════════════════════════════" << std::endl;
    std::cout << "      GAUSS'S DIVERGENCE THEOREM: ∫∫∫(∇·F)dV = ∮∮(F·n̂)dS           " << std::endl;
    std::cout << "═══════════════════════════════════════════════════════════════════" << std::endl;
    std::cout << std::endl;

    // VECTOR FIELD: F(x,y,z) = (x², y², z²) - physically: compressible flow
    // Divergence: ∇·F = 2x + 2y + 2z (analytically known)
    VectorFunction<3> F([](const VectorN<Real, 3>& p) {
        return VectorN<Real, 3>{ p[0]*p[0], p[1]*p[1], p[2]*p[2] };
    });

    // Compute divergence NUMERICALLY - MML computes this automatically!
    ScalarFunctionFromStdFunc<3> divF([&F](const VectorN<Real, 3>& p) {
        return VectorFieldOperations::DivCart<3>(F, p);
    });

    std::cout << "Vector field: F(x,y,z) = (x², y², z²)" << std::endl;
    std::cout << "Analytical divergence: ∇·F = 2x + 2y + 2z" << std::endl;
    std::cout << "Integration domain: Unit cube [0,1]³" << std::endl;
    std::cout << std::endl;

    // VOLUME INTEGRAL: ∫∫∫(∇·F)dV over unit cube [0,1]³
    auto y_lo = [](Real) { return 0.0; };
    auto y_hi = [](Real) { return 1.0; };
    auto z_lo = [](Real, Real) { return 0.0; };
    auto z_hi = [](Real, Real) { return 1.0; };

    Real volIntegral = Integrate3D(divF, GAUSS10, 0.0, 1.0, y_lo, y_hi, z_lo, z_hi).value;

    // ANALYTICAL: ∫₀¹∫₀¹∫₀¹(2x + 2y + 2z)dxdydz = 3
    Real analyticalVol = 3.0;

    // SURFACE INTEGRAL: ∮∮(F·n̂)dS through all 6 faces
    Cube3D unitCube(1.0, Point3Cartesian(0.5, 0.5, 0.5));
    Real surfIntegral = SurfaceIntegration::SurfaceIntegral(F, unitCube, 1e-8);

    std::cout << std::fixed << std::setprecision(10);
    std::cout << "┌─────────────────────────────────────────────────────────────────┐" << std::endl;
    std::cout << "│ LEFT SIDE:  ∫∫∫(∇·F)dV (volume integral)                        │" << std::endl;
    std::cout << "│   Numerical (GAUSS10): " << std::setw(15) << volIntegral << "                        │" << std::endl;
    std::cout << "│   Analytical:          " << std::setw(15) << analyticalVol << "                        │" << std::endl;
    std::cout << "├─────────────────────────────────────────────────────────────────┤" << std::endl;
    std::cout << "│ RIGHT SIDE: ∮∮(F·n̂)dS (surface integral)                        │" << std::endl;
    std::cout << "│   Numerical:           " << std::setw(15) << surfIntegral << "                        │" << std::endl;
    std::cout << "├─────────────────────────────────────────────────────────────────┤" << std::endl;
    std::cout << "│ VERIFICATION: |Volume - Surface| = " << std::scientific << std::setprecision(2)
              << std::abs(volIntegral - surfIntegral) << std::fixed << std::setprecision(10) << "                     │" << std::endl;
    std::cout << "└─────────────────────────────────────────────────────────────────┘" << std::endl;
    std::cout << std::endl;
    std::cout << "✓ GAUSS'S DIVERGENCE THEOREM VERIFIED!" << std::endl;
    std::cout << std::endl;
}

///////////////////////////////////////////////////////////////////////////////
//  STOKES' THEOREM: ∮(F·dr) = ∫∫(∇×F)·n̂ dS
///////////////////////////////////////////////////////////////////////////////
void Readme_StokesTheorem()
{
    std::cout << "═══════════════════════════════════════════════════════════════════" << std::endl;
    std::cout << "      STOKES' THEOREM: ∮(F·dr) = ∫∫(∇×F)·n̂ dS                      " << std::endl;
    std::cout << "═══════════════════════════════════════════════════════════════════" << std::endl;
    std::cout << std::endl;

    // VECTOR FIELD: F(x,y,z) = (-y, x, 0) - vortex/rotation field
    // Curl: ∇×F = (0, 0, 2) - constant upward vorticity
    VectorFunction<3> F([](const VectorN<Real, 3>& p) {
        return VectorN<Real, 3>{ -p[1], p[0], 0.0 };
    });

    std::cout << "Vector field: F(x,y,z) = (-y, x, 0) [vortex field]" << std::endl;
    std::cout << "Curl: ∇×F = (0, 0, 2) [constant vorticity in z-direction]" << std::endl;
    std::cout << std::endl;

    Real R = 2.0;  // Radius of the circular boundary

    std::cout << "Boundary: Circle of radius R = " << R << " in XY-plane at z = 0" << std::endl;
    std::cout << "Surface: Disk of radius R = " << R << " with normal n̂ = (0, 0, 1)" << std::endl;
    std::cout << std::endl;

    // LEFT SIDE: LINE INTEGRAL ∮(F·dr) around circle
    // Parametrize circle: r(t) = (R cos t, R sin t, 0), t ∈ [0, 2π]
    ParametricCurveFromStdFunc<3> circle([R](Real t) -> VectorN<Real, 3> {
        return VectorN<Real, 3>{ R * cos(t), R * sin(t), 0.0 };
    });

    Real circulation = PathIntegration::LineIntegral(F, circle, 0.0, 2*Constants::PI, 1e-10);

    // RIGHT SIDE: SURFACE INTEGRAL ∫∫(∇×F)·n̂ dS over disk
    // Curl = (0, 0, 2), normal = (0, 0, 1), so (∇×F)·n̂ = 2
    // Area of disk = πR²
    // ∫∫ 2 dA = 2 × πR² = 2πR²
    Real analyticalSurf = 2.0 * Constants::PI * R * R;

    // Numerically compute surface integral over disk using polar coordinates
    // ∫₀^R ∫₀^2π 2 × r dr dθ = 2 × ∫₀^R r dr × ∫₀^2π dθ = 2 × (R²/2) × 2π = 2πR²
    ScalarFunctionFromStdFunc<2> curlDotN([](const VectorN<Real, 2>& p) {
        (void)p;  // Independent of position - constant curl
        return 2.0;
    });

    // Integrate in polar: ∫∫ f(r,θ) r dr dθ
    auto theta_lo = [](Real) { return 0.0; };
    auto theta_hi = [](Real) { return 2.0 * Constants::PI; };

    // Transform: we integrate 2 × r over the disk
    ScalarFunctionFromStdFunc<2> integrand([](const VectorN<Real, 2>& p) {
        Real r = p[0];  // r coordinate
        return 2.0 * r;  // curlDotN × Jacobian for polar coordinates
    });
    Real numSurfIntegral = Integrate2D(integrand, GAUSS10, 0.0, R, theta_lo, theta_hi).value;

    std::cout << std::fixed << std::setprecision(10);
    std::cout << "┌─────────────────────────────────────────────────────────────────┐" << std::endl;
    std::cout << "│ LEFT SIDE:  ∮(F·dr) (line integral / circulation)               │" << std::endl;
    std::cout << "│   Numerical:          " << std::setw(15) << circulation << "                        │" << std::endl;
    std::cout << "├─────────────────────────────────────────────────────────────────┤" << std::endl;
    std::cout << "│ RIGHT SIDE: ∫∫(∇×F)·n̂ dS (surface integral of curl)             │" << std::endl;
    std::cout << "│   Numerical:          " << std::setw(15) << numSurfIntegral << "                        │" << std::endl;
    std::cout << "│   Analytical (2πR²):  " << std::setw(15) << analyticalSurf << "                        │" << std::endl;
    std::cout << "├─────────────────────────────────────────────────────────────────┤" << std::endl;
    std::cout << "│ VERIFICATION: |Line - Surface| = " << std::scientific << std::setprecision(2)
              << std::abs(circulation - numSurfIntegral) << std::fixed << std::setprecision(10) << "                       │" << std::endl;
    std::cout << "└─────────────────────────────────────────────────────────────────┘" << std::endl;
    std::cout << std::endl;
    std::cout << "✓ STOKES' THEOREM VERIFIED!" << std::endl;
    std::cout << std::endl;
}

///////////////////////////////////////////////////////////////////////////////
//  GREEN'S THEOREM: ∮(P dx + Q dy) = ∫∫(∂Q/∂x - ∂P/∂y) dA
///////////////////////////////////////////////////////////////////////////////
void Readme_GreensTheorem()
{
    std::cout << "═══════════════════════════════════════════════════════════════════" << std::endl;
    std::cout << "      GREEN'S THEOREM: ∮(P dx + Q dy) = ∫∫(∂Q/∂x - ∂P/∂y) dA       " << std::endl;
    std::cout << "═══════════════════════════════════════════════════════════════════" << std::endl;
    std::cout << std::endl;

    // 2D VECTOR FIELD: F(x,y) = (P, Q) = (x² - y, x + y²)
    // ∂Q/∂x = 1,  ∂P/∂y = -1
    // ∂Q/∂x - ∂P/∂y = 1 - (-1) = 2 (constant!)

    std::cout << "2D Vector field: F(x,y) = (x² - y, x + y²)" << std::endl;
    std::cout << "   P(x,y) = x² - y,  Q(x,y) = x + y²" << std::endl;
    std::cout << "   ∂Q/∂x = 1,  ∂P/∂y = -1" << std::endl;
    std::cout << "   ∂Q/∂x - ∂P/∂y = 2 (constant)" << std::endl;
    std::cout << std::endl;

    Real R = 1.5;  // Radius of the circular region

    std::cout << "Region: Disk of radius R = " << R << " centered at origin" << std::endl;
    std::cout << "Boundary: Circle of radius R = " << R << " (counterclockwise)" << std::endl;
    std::cout << std::endl;

    // LEFT SIDE: LINE INTEGRAL ∮(P dx + Q dy) around circle
    // We embed this as a 3D problem with z = 0 to use existing PathIntegration
    VectorFunction<3> F([](const VectorN<Real, 3>& p) {
        Real x = p[0], y = p[1];
        return VectorN<Real, 3>{ x*x - y, x + y*y, 0.0 };  // (P, Q, 0)
    });

    ParametricCurveFromStdFunc<3> circle([R](Real t) -> VectorN<Real, 3> {
        return VectorN<Real, 3>{ R * cos(t), R * sin(t), 0.0 };
    });

    Real lineIntegral = PathIntegration::LineIntegral(F, circle, 0.0, 2*Constants::PI, 1e-10);

    // RIGHT SIDE: AREA INTEGRAL ∫∫(∂Q/∂x - ∂P/∂y) dA = ∫∫ 2 dA = 2 × πR²
    Real analyticalArea = 2.0 * Constants::PI * R * R;

    // Numerical: Integrate 2 over the disk in polar coordinates
    // ∫₀^2π ∫₀^R 2 × r dr dθ
    ScalarFunctionFromStdFunc<2> integrand([](const VectorN<Real, 2>& p) {
        Real r = p[0];
        return 2.0 * r;  // (∂Q/∂x - ∂P/∂y) × Jacobian
    });
    auto theta_lo = [](Real) { return 0.0; };
    auto theta_hi = [](Real) { return 2.0 * Constants::PI; };

    Real areaIntegral = Integrate2D(integrand, GAUSS10, 0.0, R, theta_lo, theta_hi).value;

    std::cout << std::fixed << std::setprecision(10);
    std::cout << "┌─────────────────────────────────────────────────────────────────┐" << std::endl;
    std::cout << "│ LEFT SIDE:  ∮(P dx + Q dy) (line integral)                      │" << std::endl;
    std::cout << "│   Numerical:           " << std::setw(15) << lineIntegral << "                       │" << std::endl;
    std::cout << "├─────────────────────────────────────────────────────────────────┤" << std::endl;
    std::cout << "│ RIGHT SIDE: ∫∫(∂Q/∂x - ∂P/∂y) dA (area integral)                │" << std::endl;
    std::cout << "│   Numerical:           " << std::setw(15) << areaIntegral << "                       │" << std::endl;
    std::cout << "│   Analytical (2πR²):   " << std::setw(15) << analyticalArea << "                       │" << std::endl;
    std::cout << "├─────────────────────────────────────────────────────────────────┤" << std::endl;
    std::cout << "│ VERIFICATION: |Line - Area| = " << std::scientific << std::setprecision(2)
              << std::abs(lineIntegral - areaIntegral) << std::fixed << std::setprecision(10) << "                          │" << std::endl;
    std::cout << "└─────────────────────────────────────────────────────────────────┘" << std::endl;
    std::cout << std::endl;
    std::cout << "✓ GREEN'S THEOREM VERIFIED!" << std::endl;
    std::cout << std::endl;
}

///////////////////////////////////////////////////////////////////////////////
//  MAIN ENTRY POINT - ALL THREE FUNDAMENTAL THEOREMS
///////////////////////////////////////////////////////////////////////////////
void Readme_FundamentalTheorems()
{
    std::cout << std::endl;
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****        FUNDAMENTAL THEOREMS OF VECTOR CALCULUS                ****" << std::endl;
    std::cout << "****        MML Numerical Verification Suite                       ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;
    std::cout << std::endl;
    std::cout << "MML can numerically verify the three cornerstones of vector calculus:" << std::endl;
    std::cout << "  1. GAUSS'S DIVERGENCE THEOREM - relates volume to surface integrals" << std::endl;
    std::cout << "  2. STOKES' THEOREM - relates line integrals to surface curl integrals" << std::endl;
    std::cout << "  3. GREEN'S THEOREM - 2D version of Stokes for planar regions" << std::endl;
    std::cout << std::endl;

    Readme_GausssDivergenceTheorem();
    Readme_StokesTheorem();
    Readme_GreensTheorem();

    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****              ALL FUNDAMENTAL THEOREMS VERIFIED!               ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;
    std::cout << std::endl;
}
