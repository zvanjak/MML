///////////////////////////////////////////////////////////////////////////////////////////
/// @file parametric_curves_visualization.cpp
/// @brief Chapter 02: Parametric Curve Visualization Demos
/// @details Demonstrates visualization of parametric curves including:
///          - 2D curves: Butterfly curve and other beautiful mathematical curves
///          - 3D curves: Lorenz attractor, space curves, and helices
///////////////////////////////////////////////////////////////////////////////////////////

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "mml/base/Vector/VectorN.h"
#include "mml/base/Function.h"
#include "mml/base/InterpolatedFunction.h"
#include "mml/core/FunctionHelpers.h"
#include "mml/core/Curves.h"
#include "mml/base/ODESystem.h"
#include "mml/algorithms/ODESolvers/ODEAdaptiveIntegrator.h"
#include "mml/tools/Visualizer.h"
#endif

using namespace MML;
using namespace MML::Curves;

///////////////////////////////////////////////////////////////////////////////////////////
/// @brief Demo: 2D Parametric Curves - The Butterfly Curve
/// @details The Butterfly curve is a beautiful transcendental plane curve discovered
///          by Temple H. Fay. The parametric equations are:
///          x(t) = sin(t) * (e^cos(t) - 2*cos(4t) - sin^5(t/12))
///          y(t) = cos(t) * (e^cos(t) - 2*cos(4t) - sin^5(t/12))
///          
///          Also demonstrates other classic 2D parametric curves from the library.
///////////////////////////////////////////////////////////////////////////////////////////
void Chapter02_Demo_ParametricCurves2D()
{
    std::cout << "\n=== Demo: 2D Parametric Curves ===\n";
    
    //-------------------------------------------------------------------------
    // The Butterfly Curve - Temple H. Fay's beautiful discovery
    // Now using the predefined ButterflyCurve class from Curves.h
    //-------------------------------------------------------------------------
    ButterflyCurve butterfly;
    
    Visualizer::VisualizeParamCurve2D(
        butterfly, "Butterfly Curve (Temple H. Fay)",
        0.0, 12.0 * Constants::PI, 2001, "ch02_butterfly_curve.mml");
    
    std::cout << "   Visualized: Butterfly curve - a transcendental beauty!\n";
    
    //-------------------------------------------------------------------------
    // Classic library curves for comparison
    //-------------------------------------------------------------------------
    
    // Lemniscate of Bernoulli - the infinity symbol
    LemniscateCurve lemniscate;
    Visualizer::VisualizeParamCurve2D(
        lemniscate, "Lemniscate of Bernoulli (∞)",
        0.0, 2.0 * Constants::PI, 201, "ch02_lemniscate.mml");
    
    // Deltoid - the three-cusped hypocycloid
    DeltoidCurve deltoid(5.0);
    Visualizer::VisualizeParamCurve2D(
        deltoid, "Deltoid (3-cusped hypocycloid)",
        0.0, 2.0 * Constants::PI, 201, "ch02_deltoid.mml");
    
    // Astroid - the four-cusped hypocycloid
    AstroidCurve astroid(5.0);
    Visualizer::VisualizeParamCurve2D(
        astroid, "Astroid (4-cusped hypocycloid)",
        0.0, 2.0 * Constants::PI, 201, "ch02_astroid.mml");
    
    // Logarithmic spiral family - nature's favorite curve
    LogSpiralCurve spiral1(-0.15), spiral2(-0.25), spiral3(-0.35);
    Visualizer::VisualizeMultiParamCurve2D(
        { &spiral1, &spiral2, &spiral3 },
        "Logarithmic Spirals (Nature's Curve)",
        0.0, 20.0, 501, "ch02_log_spirals.mml");
    
    std::cout << "   Also visualized: Lemniscate, Deltoid, Astroid, Log Spirals\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
/// @brief Demo: 3D Parametric Curves - Lorenz Attractor
/// @details The famous Lorenz system - a classic example of deterministic chaos.
///          Solved as an ODE system, then visualized as a 3D parametric curve.
///          σ=10, ρ=28, β=8/3 (classic parameters)
///          
///          This demonstrates the "butterfly effect" - sensitive dependence on ICs.
///////////////////////////////////////////////////////////////////////////////////////////
void Chapter02_Demo_ParametricCurves3D()
{
    std::cout << "\n=== Demo: 3D Parametric Curves - Lorenz Attractor ===\n";
    
    //-------------------------------------------------------------------------
    // The Lorenz System - discovered by Edward Lorenz in 1963
    // Models atmospheric convection, exhibits chaotic behavior
    //-------------------------------------------------------------------------
    ODESystem lorenz_system(3, [](Real t, const Vector<Real>& x, Vector<Real>& dxdt)
    {
        const Real sigma = 10.0;   // Prandtl number
        const Real rho = 28.0;     // Rayleigh number
        const Real beta = 8.0 / 3.0;
        
        dxdt[0] = sigma * (x[1] - x[0]);           // dx/dt
        dxdt[1] = x[0] * (rho - x[2]) - x[1];      // dy/dt  
        dxdt[2] = x[0] * x[1] - beta * x[2];       // dz/dt
    });
    
    Vector<Real> initial_state({ 1.0, 1.0, 1.0 });
    
    CashKarpIntegrator solver(lorenz_system);
    ODESystemSolution sol = solver.integrate(initial_state, 0.0, 50.0, 0.001, 1e-10, 0.001);
    
    // Scale factor for visualization - makes the attractor bigger on screen
    const Real scale = 5.0;
    
    // Extract solution components and scale them
    Vector<Real> t_vals = sol.getTValues();
    Vector<Real> x_scaled = sol.getXValues(0) * scale;
    Vector<Real> y_scaled = sol.getXValues(1) * scale;
    Vector<Real> z_scaled = sol.getXValues(2) * scale;
    
    // Create scaled parametric curve using spline interpolation
    SplineInterpRealFunc x_spline(t_vals, x_scaled);
    SplineInterpRealFunc y_spline(t_vals, y_scaled);
    SplineInterpRealFunc z_spline(t_vals, z_scaled);
    
    ParametricCurveFromStdFunc<3> scaled_attractor(0.0, 50.0, [&](Real t) {
        return VectorN<Real, 3>({ x_spline(t), y_spline(t), z_spline(t) });
    });
    
    // Visualize as 3D parametric curve - the famous butterfly attractor
    Visualizer::VisualizeParamCurve3D(
        scaled_attractor, "Lorenz Attractor (3D) [scaled x5]", 
        0.0, 50.0, 5000, "ch02_lorenz_attractor_3d.mml");
    
    std::cout << "   Visualized: Lorenz attractor as 3D parametric curve\n";
    
    //-------------------------------------------------------------------------
    // Butterfly Effect - three trajectories from slightly different ICs
    //-------------------------------------------------------------------------
    Vector<Real> ic1({ 1.0, 1.0, 1.0 });
    Vector<Real> ic2({ 1.001, 1.0, 1.0 });  // Tiny perturbation in x!
    Vector<Real> ic3({ 1.0, 1.001, 1.0 });  // Tiny perturbation in y!
    
    ODESystemSolution sol1 = solver.integrate(ic1, 0.0, 30.0, 0.001, 1e-10, 0.001);
    ODESystemSolution sol2 = solver.integrate(ic2, 0.0, 30.0, 0.001, 1e-10, 0.001);
    ODESystemSolution sol3 = solver.integrate(ic3, 0.0, 30.0, 0.001, 1e-10, 0.001);
    
    // Scale all three solutions - use sx/sy/sz naming to avoid conflict with Bessel functions
    Vector<Real> tv1 = sol1.getTValues(), tv2 = sol2.getTValues(), tv3 = sol3.getTValues();
    SplineInterpRealFunc sx1(tv1, sol1.getXValues(0) * scale), sy1(tv1, sol1.getXValues(1) * scale), sz1(tv1, sol1.getXValues(2) * scale);
    SplineInterpRealFunc sx2(tv2, sol2.getXValues(0) * scale), sy2(tv2, sol2.getXValues(1) * scale), sz2(tv2, sol2.getXValues(2) * scale);
    SplineInterpRealFunc sx3(tv3, sol3.getXValues(0) * scale), sy3(tv3, sol3.getXValues(1) * scale), sz3(tv3, sol3.getXValues(2) * scale);
    
    ParametricCurveFromStdFunc<3> curve1(0.0, 30.0, [&](Real t) { return VectorN<Real, 3>({ sx1(t), sy1(t), sz1(t) }); });
    ParametricCurveFromStdFunc<3> curve2(0.0, 30.0, [&](Real t) { return VectorN<Real, 3>({ sx2(t), sy2(t), sz2(t) }); });
    ParametricCurveFromStdFunc<3> curve3(0.0, 30.0, [&](Real t) { return VectorN<Real, 3>({ sx3(t), sy3(t), sz3(t) }); });
    
    Visualizer::VisualizeMultiParamCurve3D(
        { &curve1, &curve2, &curve3 }, 
        "Lorenz: Butterfly Effect - Chaos in Action [scaled x5]", 
        0.0, 30.0, 5000, "ch02_lorenz_butterfly_effect_3d.mml");
    
    std::cout << "   The 'butterfly effect': tiny changes (0.001) → vastly different paths!\n";
    
    //-------------------------------------------------------------------------
    // Bonus: Classic 3D curves from the library
    //-------------------------------------------------------------------------
    
    // Helix - the simplest 3D curve
    HelixCurve helix(5.0, 2.0);  // radius=5, pitch parameter b=2
    Visualizer::VisualizeParamCurve3D(
        helix, "Helix (DNA, springs, screws)",
        0.0, 6.0 * Constants::PI, 501, "ch02_helix.mml");
    
    // Toroidal spiral - curve winding around a torus
    ToroidalSpiralCurve toroidal(5, 3.0);  // 5 windings, scale=3
    Visualizer::VisualizeParamCurve3D(
        toroidal, "Toroidal Spiral (5 windings)",
        0.0, 2.0 * Constants::PI, 501, "ch02_toroidal_spiral.mml");
    
    std::cout << "   Also visualized: Helix, Toroidal Spiral\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
/// @brief Entry point for Parametric Curve visualization demos
///////////////////////////////////////////////////////////////////////////////////////////
void Chapter02_ParametricCurveVisualization()
{
    std::cout << "\n";
    std::cout << "╔═══════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║         CHAPTER 02: Parametric Curve Visualization                ║\n";
    std::cout << "╚═══════════════════════════════════════════════════════════════════╝\n";
    
    Chapter02_Demo_ParametricCurves2D();
    Chapter02_Demo_ParametricCurves3D();
}
