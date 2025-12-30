/**
 * @file docs_demo_visualizers.cpp
 * @brief Demonstration code for Visualizer documentation examples
 * 
 * This file contains working implementations of all examples from docs/tools/Visualizers.md
 * Each function demonstrates a specific visualization capability and verifies the documentation
 * examples actually compile and work correctly.
 */

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/Function.h"
#include "base/VectorN.h"
#include "core/Fields.h"
#include "core/Curves.h"
#include "core/FunctionHelpers.h"

#include "base/ODESystem.h"
#include "base/ODESystemSolution.h"
#include "algorithms/ODESystemSolver.h"
#include "algorithms/ODESystemStepCalculators.h"
#include "algorithms/ODESystemSteppers.h"

#include "tools/Visualizer.h"
#include "tools/Serializer.h"
#endif

using namespace MML;

///////////////////////////////////////////////////////////////////////////////////////////
///                    EXAMPLE 1: FUNCTION ANALYSIS                                     ///
///////////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Visualize a Gaussian function and compare with Taylor approximations
 * 
 * From Visualizers.md - Example 1: Function Analysis
 */
void Docs_Demo_Visualizers_Example1_FunctionAnalysis()
{
    std::cout << "\n==========================================================================\n";
    std::cout << "Example 1: Function Analysis - Gaussian and Taylor Approximations\n";
    std::cout << "==========================================================================\n";

    // Define Gaussian function
    auto gaussian = [](Real x) {
        return std::exp(-x*x);
    };
    RealFunction f(gaussian);
    
    std::cout << "1a. Visualizing Gaussian exp(-x^2) over [-3, 3]\n";
    auto result1 = Visualizer::VisualizeRealFunction(
        f,
        "Gaussian exp(-x^2)",
        -3.0, 3.0,
        200,
        "docs_gaussian.txt"
    );
    std::cout << "    Result: " << (result1.success ? "Success" : "Failed: " + result1.errorMessage) << "\n";
    
    // Taylor approximations of exp(-x^2) around x=0:
    // O(2): 1 - x^2
    // O(4): 1 - x^2 + x^4/2
    auto taylor2 = [](Real x) { return 1.0 - x*x; };
    auto taylor4 = [](Real x) { return 1.0 - x*x + x*x*x*x/2.0; };
    
    RealFunction t2(taylor2);
    RealFunction t4(taylor4);
    std::vector<IRealFunction*> funcs = {&f, &t2, &t4};
    std::vector<std::string> legend = {"exp(-x^2)", "Taylor O(2)", "Taylor O(4)"};
    
    std::cout << "1b. Comparing Gaussian with Taylor approximations\n";
    auto result2 = Visualizer::VisualizeMultiRealFunction(
        funcs,
        "Taylor Approximations of Gaussian",
        legend,
        -2.0, 2.0,
        150,
        "docs_taylor_comparison.txt"
    );
    std::cout << "    Result: " << (result2.success ? "Success" : "Failed: " + result2.errorMessage) << "\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                    EXAMPLE 2: VECTOR FIELD VISUALIZATION                            ///
///////////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Visualize electric field of a point charge
 * 
 * From Visualizers.md - Example 2: Vector Field Visualization
 */
void Docs_Demo_Visualizers_Example2_VectorField()
{
    std::cout << "\n==========================================================================\n";
    std::cout << "Example 2: Vector Field - Electric Field of Point Charge\n";
    std::cout << "==========================================================================\n";

    // Electric field of point charge at origin: E = r / |r|^3
    auto eField = [](const VectorN<Real, 2>& r) {
        Real x = r[0], y = r[1];
        Real rMag = std::sqrt(x*x + y*y);
        if (rMag < 1e-10) return VectorN<Real, 2>{0.0, 0.0};
        
        Real scale = 1.0 / (rMag * rMag * rMag);
        return VectorN<Real, 2>{x * scale, y * scale};
    };
    
    VectorFunction<2> field(eField);
    
    std::cout << "Visualizing electric field of point charge at origin\n";
    std::cout << "Field: E = r / |r|^3 (inverse-square law)\n";
    auto result = Visualizer::VisualizeVectorField2DCartesian(
        field,
        "Electric Field (Point Charge)",
        -5.0, 5.0, 25,
        -5.0, 5.0, 25,
        "docs_electric_field.txt"
    );
    std::cout << "Result: " << (result.success ? "Success" : "Failed: " + result.errorMessage) << "\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                    EXAMPLE 3: ODE SOLUTION PHASE PORTRAIT                           ///
///////////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Visualize simple pendulum: time series, phase portrait, and multi-function
 * 
 * From Visualizers.md - Example 3: ODE Solution Phase Portrait
 */
void Docs_Demo_Visualizers_Example3_ODEPendulum()
{
    std::cout << "\n==========================================================================\n";
    std::cout << "Example 3: ODE Solution - Simple Pendulum Phase Portrait\n";
    std::cout << "==========================================================================\n";

    // Pendulum: θ'' + sin(θ) = 0
    // State: [θ, ω] where ω = θ'
    // Equations: dθ/dt = ω, dω/dt = -sin(θ)
    auto pendulumFunc = [](Real t, const Vector<Real>& state, Vector<Real>& dstate) {
        Real theta = state[0];
        Real omega = state[1];
        dstate[0] = omega;
        dstate[1] = -std::sin(theta);
    };
    
    ODESystem pendulum(2, pendulumFunc);
    Vector<Real> initState({Constants::PI/4.0, 0.0});  // 45° initial angle, zero velocity
    
    std::cout << "System: Simple pendulum θ'' + sin(θ) = 0\n";
    std::cout << "Initial conditions: θ = 45°, ω = 0\n";
    std::cout << "Integration: t ∈ [0, 10] with 1000 steps (RK4)\n\n";
    
    // Solve using RK4 fixed step
    ODESystemFixedStepSolver solver(pendulum, StepCalculators::RK4_Basic);
    ODESystemSolution sol = solver.integrate(initState, 0.0, 10.0, 1000);
    
    std::cout << "Solution obtained: " << sol.getTotalSavedSteps() << " steps saved\n\n";

    // 3a. Time series of angle
    std::cout << "3a. Visualizing angle θ(t) vs time\n";
    auto result1 = Visualizer::VisualizeODESysSolCompAsFunc(
        sol, 0, "Pendulum Angle vs Time", "docs_theta_t.txt"
    );
    std::cout << "    Result: " << (result1.success ? "Success" : "Failed: " + result1.errorMessage) << "\n";
    
    // 3b. Phase portrait (ω vs θ)
    std::cout << "3b. Visualizing phase portrait (ω vs θ)\n";
    auto result2 = Visualizer::VisualizeODESysSolAsParamCurve2(
        sol, 0, 1, "Pendulum Phase Portrait", "docs_phase.txt"
    );
    std::cout << "    Result: " << (result2.success ? "Success" : "Failed: " + result2.errorMessage) << "\n";
    
    // 3c. Both variables as multi-function
    std::cout << "3c. Visualizing both variables θ(t) and ω(t)\n";
    std::vector<std::string> labels = {"θ (angle)", "ω (angular velocity)"};
    auto result3 = Visualizer::VisualizeODESysSolAsMultiFunc(
        sol, "Pendulum State Variables", labels, "docs_pendulum_state.txt"
    );
    std::cout << "    Result: " << (result3.success ? "Success" : "Failed: " + result3.errorMessage) << "\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                    EXAMPLE 4: 3D PARAMETRIC CURVES                                  ///
///////////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Visualize trefoil knot as a 3D parametric curve
 * 
 * From Visualizers.md - Example 4: 3D Curves
 */
void Docs_Demo_Visualizers_Example4_3DCurves()
{
    std::cout << "\n==========================================================================\n";
    std::cout << "Example 4: 3D Parametric Curves - Trefoil Knot\n";
    std::cout << "==========================================================================\n";

    // Trefoil knot parametrization:
    // x(t) = sin(t) + 2*sin(2t)
    // y(t) = cos(t) - 2*cos(2t)
    // z(t) = -sin(3t)
    // Scaled by 50 for proper visualization display
    auto trefoil = [](Real t) {
        Real x = 50.0 * (std::sin(t) + 2.0*std::sin(2.0*t));
        Real y = 50.0 * (std::cos(t) - 2.0*std::cos(2.0*t));
        Real z = 50.0 * (-std::sin(3.0*t));
        return VectorN<Real, 3>{x, y, z};
    };
    
    ParametricCurve<3> knot(trefoil);
    
    std::cout << "Trefoil Knot:\n";
    std::cout << "  x(t) = sin(t) + 2*sin(2t)\n";
    std::cout << "  y(t) = cos(t) - 2*cos(2t)\n";
    std::cout << "  z(t) = -sin(3t)\n";
    std::cout << "  t ∈ [0, 2π], 500 points\n\n";
    
    auto result = Visualizer::VisualizeParamCurve3D(
        knot,
        "Trefoil Knot",
        0.0, 2.0*Constants::PI,
        500,
        "docs_trefoil.txt"
    );
    std::cout << "Result: " << (result.success ? "Success" : "Failed: " + result.errorMessage) << "\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                    ADDITIONAL EXAMPLES                                              ///
///////////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Visualize scalar function as 3D surface
 */
void Docs_Demo_Visualizers_ScalarFunction()
{
    std::cout << "\n==========================================================================\n";
    std::cout << "Additional: Scalar Function as 3D Surface\n";
    std::cout << "==========================================================================\n";

    // f(x,y) = sin(x) * cos(y)
    // Scaled by 25 for proper visualization display
    auto surfFunc = [](const VectorN<Real, 2>& v) {
        return 25.0 * std::sin(v[0]) * std::cos(v[1]);
    };
    ScalarFunction<2> surface(surfFunc);

    std::cout << "Visualizing f(x,y) = sin(x)*cos(y)\n";
    std::cout << "Domain: x,y ∈ [-π, π], 50x50 grid\n\n";
    
    auto result = Visualizer::VisualizeScalarFunc2DCartesian(
        surface,
        "sin(x)*cos(y)",
        -Constants::PI, Constants::PI, 50,
        -Constants::PI, Constants::PI, 50,
        "docs_wave_surface.txt"
    );
    std::cout << "Result: " << (result.success ? "Success" : "Failed: " + result.errorMessage) << "\n";
}

/**
 * @brief Visualize 2D parametric curve (Lissajous)
 */
void Docs_Demo_Visualizers_ParamCurve2D()
{
    std::cout << "\n==========================================================================\n";
    std::cout << "Additional: 2D Parametric Curve - Lissajous Figure\n";
    std::cout << "==========================================================================\n";

    // Lissajous curve: x = sin(3t), y = sin(4t)
    auto lissajous = [](Real t) {
        return VectorN<Real, 2>{std::sin(3.0*t), std::sin(4.0*t)};
    };
    ParametricCurve<2> curve(lissajous);

    std::cout << "Lissajous curve: x = sin(3t), y = sin(4t)\n";
    std::cout << "t ∈ [0, 2π], 500 points\n\n";
    
    auto result = Visualizer::VisualizeParamCurve2D(
        curve,
        "Lissajous 3:4",
        0.0, 2.0*Constants::PI,
        500,
        "docs_lissajous.txt"
    );
    std::cout << "Result: " << (result.success ? "Success" : "Failed: " + result.errorMessage) << "\n";
}

/**
 * @brief Visualize multiple functions comparison
 */
void Docs_Demo_Visualizers_MultiFunctions()
{
    std::cout << "\n==========================================================================\n";
    std::cout << "Additional: Multiple Function Comparison\n";
    std::cout << "==========================================================================\n";

    // Compare sin, cos, and their product
    RealFunction sine([](Real x) { return std::sin(x); });
    RealFunction cosine([](Real x) { return std::cos(x); });
    RealFunction product([](Real x) { return std::sin(x) * std::cos(x); });

    std::vector<IRealFunction*> funcs = {&sine, &cosine, &product};
    std::vector<std::string> legend = {"sin(x)", "cos(x)", "sin(x)*cos(x)"};

    std::cout << "Comparing: sin(x), cos(x), sin(x)*cos(x)\n";
    std::cout << "x ∈ [0, 2π], 200 points\n\n";
    
    auto result = Visualizer::VisualizeMultiRealFunction(
        funcs,
        "Trigonometric Function Comparison",
        legend,
        0.0, 2.0*Constants::PI,
        200,
        "docs_trig_compare.txt"
    );
    std::cout << "Result: " << (result.success ? "Success" : "Failed: " + result.errorMessage) << "\n";
}

/**
 * @brief Visualize Lorenz attractor as 3D parametric curve
 */
void Docs_Demo_Visualizers_LorenzAttractor()
{
    std::cout << "\n==========================================================================\n";
    std::cout << "Additional: Lorenz Attractor as 3D Curve\n";
    std::cout << "==========================================================================\n";

    // Lorenz system: dx/dt = σ(y-x), dy/dt = x(ρ-z)-y, dz/dt = xy-βz
    // Classic parameters: σ=10, ρ=28, β=8/3
    // Using constants since ODESystem requires a function pointer (no captures)
    // State scaled by 5 for proper visualization display
    auto lorenzFunc = [](Real t, const Vector<Real>& state, Vector<Real>& dstate) {
        constexpr Real sigma = 10.0;
        constexpr Real rho = 28.0;
        constexpr Real beta = 8.0/3.0;
        constexpr Real scale = 5.0;
        Real x = state[0]/scale, y = state[1]/scale, z = state[2]/scale;
        dstate[0] = scale * sigma * (y - x);
        dstate[1] = scale * (x * (rho - z) - y);
        dstate[2] = scale * (x * y - beta * z);
    };
    
    ODESystem lorenz(3, lorenzFunc);
    Vector<Real> initState({1.0, 1.0, 1.0});
    
    std::cout << "Lorenz system: σ=10, ρ=28, β=8/3\n";
    std::cout << "Initial state: (1, 1, 1)\n";
    std::cout << "Integration: t ∈ [0, 50] with 10000 steps (RK4)\n\n";
    
    ODESystemFixedStepSolver solver(lorenz, StepCalculators::RK4_Basic);
    ODESystemSolution sol = solver.integrate(initState, 0.0, 50.0, 10000);
    
    std::cout << "Solution obtained: " << sol.getTotalSavedSteps() << " steps\n\n";

    auto result = Visualizer::VisualizeODESysSolAsParamCurve3(
        sol, 0, 1, 2, "Lorenz Attractor", "docs_lorenz_attractor.txt"
    );
    std::cout << "Result: " << (result.success ? "Success" : "Failed: " + result.errorMessage) << "\n";
}

/**
 * @brief Visualize 3D vector field
 */
void Docs_Demo_Visualizers_VectorField3D()
{
    std::cout << "\n==========================================================================\n";
    std::cout << "Additional: 3D Vector Field\n";
    std::cout << "==========================================================================\n";

    // Simple 3D field: F = (y, z, x) - cyclic rotation
    auto field3D = [](const VectorN<Real, 3>& r) {
        return VectorN<Real, 3>{r[1], r[2], r[0]};
    };
    VectorFunction<3> field(field3D);

    std::cout << "3D Vector Field: F(x,y,z) = (y, z, x)\n";
    std::cout << "Domain: [-3,3]^3, 8x8x8 grid\n\n";
    
    auto result = Visualizer::VisualizeVectorField3DCartesian(
        field,
        "Cyclic Rotation Field",
        -3.0, 3.0, 8,
        -3.0, 3.0, 8,
        -3.0, 3.0, 8,
        "docs_vector3d_cyclic.txt"
    );
    std::cout << "Result: " << (result.success ? "Success" : "Failed: " + result.errorMessage) << "\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                    MASTER DEMO FUNCTION                                             ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_Visualizers()
{
    std::cout << "***********************************************************" << std::endl;
    std::cout << "*****           VISUALIZERS DOCUMENTATION DEMOS       *****" << std::endl;
    std::cout << "***********************************************************" << std::endl;
    std::cout << "\nThese examples verify that all code in Visualizers.md works correctly.\n";
    std::cout << "Output files are created in the 'results/' directory.\n";

    // Main examples from documentation
    Docs_Demo_Visualizers_Example1_FunctionAnalysis();
    Docs_Demo_Visualizers_Example2_VectorField();
    Docs_Demo_Visualizers_Example3_ODEPendulum();
    Docs_Demo_Visualizers_Example4_3DCurves();
    
    // Additional examples
    Docs_Demo_Visualizers_ScalarFunction();
    Docs_Demo_Visualizers_ParamCurve2D();
    Docs_Demo_Visualizers_MultiFunctions();
    Docs_Demo_Visualizers_LorenzAttractor();
    Docs_Demo_Visualizers_VectorField3D();

    std::cout << "\n***********************************************************" << std::endl;
    std::cout << "*****           ALL VISUALIZER DEMOS COMPLETE          *****" << std::endl;
    std::cout << "***********************************************************" << std::endl;
    std::cout << "\nGenerated files in results/:\n";
    std::cout << "  - docs_gaussian.txt\n";
    std::cout << "  - docs_taylor_comparison.txt\n";
    std::cout << "  - docs_electric_field.txt\n";
    std::cout << "  - docs_theta_t.txt\n";
    std::cout << "  - docs_phase.txt\n";
    std::cout << "  - docs_pendulum_state.txt\n";
    std::cout << "  - docs_trefoil.txt\n";
    std::cout << "  - docs_wave_surface.txt\n";
    std::cout << "  - docs_lissajous.txt\n";
    std::cout << "  - docs_trig_compare.txt\n";
    std::cout << "  - docs_lorenz_attractor.txt\n";
    std::cout << "  - docs_vector3d_cyclic.txt\n";
}
