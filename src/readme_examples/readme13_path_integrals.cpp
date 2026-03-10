///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        readme13_path_integrals.cpp                                         ///
///  Description: README example - Path and Line Integrals                            ///
///               Demonstrates work integral, arc length, path integration            ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///////////////////////////////////////////////////////////////////////////////////////////

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/Vector/VectorN.h"
#include "base/Function.h"
#include "core/Integration/PathIntegration.h"
#endif

#include <iostream>
#include <iomanip>

using namespace MML;

void Readme_PathIntegrals()
{
    std::cout << std::endl;
    std::cout << "=== Path & Line Integrals ===" << std::endl;

    // Line integral of vector field along a curve
    // ∫ F·dr where F = (y, -x, 0) along unit circle
    // This is a rotational field - work around closed loop is non-zero!

    VectorFunction<3> F([](const VectorN<Real, 3>& p) -> VectorN<Real, 3> {
        return {p[1], -p[0], 0.0};
    });

    // Unit circle in XY plane: r(t) = (cos(t), sin(t), 0), t ∈ [0, 2π]
    ParametricCurve<3> circle([](Real t) -> VectorN<Real, 3> {
        return {cos(t), sin(t), 0.0};
    });

    std::cout << "Vector field F = (y, -x, 0)" << std::endl;
    std::cout << "Path: unit circle in XY plane" << std::endl;

    // Work integral: W = ∫₀^{2π} F(r(t)) · r'(t) dt
    Real work = PathIntegration::LineIntegral(F, circle, 0.0, 2*Constants::PI, 1e-8);
    std::cout << std::endl << "Work integral ∮ F·dr:" << std::endl;
    std::cout << "  Numerical:  " << work << std::endl;
    std::cout << "  Analytical: -2π = " << -2*Constants::PI << std::endl;

    // Arc length computation: ∫ ds
    Real arc_length = PathIntegration::ParametricCurveLength<3>(circle, 0.0, 2*Constants::PI);
    std::cout << std::endl << "Arc length of unit circle:" << std::endl;
    std::cout << "  Numerical:  " << arc_length << std::endl;
    std::cout << "  Analytical: 2π = " << 2*Constants::PI << std::endl;

    // Scalar field line integral: ∫ f ds (arc length weighting)
    // Integrate temperature field T(x,y,z) = x² + y² along circle
    ScalarFunction<3> temperature([](const VectorN<Real, 3>& p) { 
        return p[0]*p[0] + p[1]*p[1];  // r² = 1 on unit circle
    });
    
    Real temp_integral = PathIntegration::LineIntegral(temperature, circle, 
                                                        0.0, 2*Constants::PI, 1e-8);
    std::cout << std::endl << "Scalar line integral ∫ (x² + y²) ds on unit circle:" << std::endl;
    std::cout << "  Numerical:  " << temp_integral << std::endl;
    std::cout << "  Expected:   " << 2*Constants::PI << " (since x² + y² = 1 everywhere)" << std::endl;

    // Semi-circle path
    Real half_work = PathIntegration::LineIntegral(F, circle, 0.0, Constants::PI, 1e-8);
    std::cout << std::endl << "Work over half circle [0,π]: " << half_work << std::endl;
    std::cout << "  Expected: -π = " << -Constants::PI << std::endl;

    std::cout << std::endl;
}
