///////////////////////////////////////////////////////////////////////////////////////////
/// @file chapter02_real_functions.cpp
/// @brief Chapter 02: Real Function Visualization Demos
/// @details Demonstrates visualization of real-valued functions including:
///          - Classic mathematical functions (Sinc, Gaussian, Damped oscillation)
///          - Functions with their derivatives
///          - Lorenz attractor and chaos visualization
///////////////////////////////////////////////////////////////////////////////////////////

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "mml/base/Vector/VectorN.h"
#include "mml/base/Function.h"
#include "mml/base/InterpolatedFunction.h"
#include "mml/core/FunctionHelpers.h"
#include "mml/base/ODESystem.h"
#include "mml/algorithms/ODESolvers/ODEAdaptiveIntegrator.h"
#include "mml/tools/Visualizer.h"
#endif

using namespace MML;

///////////////////////////////////////////////////////////////////////////////////////////
/// @brief Demo 1: Classic Mathematical Functions Gallery
/// @details Shows three fundamental functions together:
///          - Sinc function: sin(x)/x - fundamental in signal processing, Fourier analysis
///          - Gaussian: exp(-x²/2) - probability, quantum mechanics, heat diffusion
///          - Damped oscillation: exp(-|x|/3)*cos(x) - physics, engineering, decay
///////////////////////////////////////////////////////////////////////////////////////////
void Chapter02_Demo_ClassicFunctions()
{
    std::cout << "\n=== Demo 1: Classic Mathematical Functions ===\n";
    
    // Sinc function - cardinal sine, fundamental in signal processing
    // Note: handle x=0 specially to avoid division by zero
    RealFunctionFromStdFunc sinc{ [](Real x) { 
        return std::abs(x) < 1e-10 ? Real(1.0) : std::sin(x) / x; 
    }};
    
    // Gaussian / Bell curve - appears everywhere in probability and physics
    RealFunctionFromStdFunc gaussian{ [](Real x) { 
        return std::exp(-x * x / 2.0); 
    }};
    
    // Damped oscillation - models decay in physical systems
    RealFunctionFromStdFunc damped{ [](Real x) { 
        return std::exp(-std::abs(x) / 3.0) * std::cos(x); 
    }};
    
    // Visualize each individually
    Visualizer::VisualizeRealFunction(sinc, "Sinc Function: sin(x)/x", 
                                      -15.0, 15.0, 501, "ch02_sinc.mml");
    
    Visualizer::VisualizeRealFunction(gaussian, "Gaussian: exp(-x²/2)", 
                                      -5.0, 5.0, 301, "ch02_gaussian.mml");
    
    Visualizer::VisualizeRealFunction(damped, "Damped Oscillation: exp(-|x|/3)·cos(x)", 
                                      -15.0, 15.0, 501, "ch02_damped.mml");
    
    // Show all three together for comparison
    Visualizer::VisualizeMultiRealFunction(
        { &sinc, &gaussian, &damped },
        "Classic Functions Gallery",
        { "Sinc: sin(x)/x", "Gaussian: exp(-x²/2)", "Damped: exp(-|x|/3)·cos(x)" },
        -10.0, 10.0, 501, "ch02_classic_functions.mml");
    
    std::cout << "   Visualized: Sinc, Gaussian, Damped Oscillation\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
/// @brief Demo 2: Function with Derivatives
/// @details Shows a function together with its first and second derivatives.
///          Using f(x) = sin(x)·exp(-x²/10) - a modulated Gaussian
///          f'(x) and f''(x) computed numerically using MML's differentiation
///////////////////////////////////////////////////////////////////////////////////////////
void Chapter02_Demo_FunctionWithDerivatives()
{
    std::cout << "\n=== Demo 2: Function with First and Second Derivatives ===\n";
    
    // A nice function: modulated Gaussian - has clear features
    RealFunctionFromStdFunc f{ [](Real x) { 
        return std::sin(x) * std::exp(-x * x / 10.0); 
    }};
    
    // Create numerical derivatives using MML's derivation classes
    // RealFuncDerived4 = 1st derivative with 4th order accuracy
    // RealFuncSecondDerived4 = 2nd derivative with 4th order accuracy
    RealFuncDerived4 f_prime(f);              // First derivative
    RealFuncSecondDerived4 f_double_prime(f); // Second derivative
    
    // Visualize function alone
    Visualizer::VisualizeRealFunction(f, "f(x) = sin(x)·exp(-x²/10)", 
                                      -8.0, 8.0, 401, "ch02_func.mml");
    
    // Visualize first derivative
    Visualizer::VisualizeRealFunction(f_prime, "f'(x) - First Derivative", 
                                      -8.0, 8.0, 401, "ch02_func_deriv1.mml");
    
    // Visualize second derivative
    Visualizer::VisualizeRealFunction(f_double_prime, "f''(x) - Second Derivative", 
                                      -8.0, 8.0, 401, "ch02_func_deriv2.mml");
    
    // Show all three together - the educational visualization
    Visualizer::VisualizeMultiRealFunction(
        { &f, &f_prime, &f_double_prime },
        "Function and Its Derivatives",
        { "f(x) = sin(x)·exp(-x²/10)", "f'(x) - First Derivative", "f''(x) - Second Derivative" },
        -8.0, 8.0, 401, "ch02_func_with_derivatives.mml");
    
    std::cout << "   Visualized: f(x), f'(x), f''(x)\n";
    std::cout << "   Notice: f'(x)=0 at extrema, f''(x)=0 at inflection points\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
/// @brief Demo 3: Lorenz System Time Series
/// @details The famous Lorenz system - visualized as real functions x(t), y(t), z(t).
///          Shows the chaotic time evolution of the three state variables.
///          σ=10, ρ=28, β=8/3 (classic parameters)
///          
///          Note: The 3D attractor visualization is in parametric_curves_visualization.cpp
///////////////////////////////////////////////////////////////////////////////////////////
void Chapter02_Demo_LorenzTimeSeries()
{
    std::cout << "\n=== Demo 3: Lorenz System as Time Series ===\n";
    
    // The Lorenz system - discovered by Edward Lorenz in 1963
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
    
    // Extract solution components as real functions of time
    Vector<Real> t_vals = sol.getTValues();
    SplineInterpRealFunc x_t(t_vals, sol.getXValues(0));
    SplineInterpRealFunc y_t(t_vals, sol.getXValues(1));
    SplineInterpRealFunc z_t(t_vals, sol.getXValues(2));
    
    // Visualize as time series - each component as a real function
    std::vector<IRealFunction*> time_series = { &x_t, &y_t, &z_t };
    Visualizer::VisualizeMultiRealFunction(
        time_series,
        "Lorenz System: x(t), y(t), z(t) - Chaotic Time Evolution",
        { "x(t)", "y(t)", "z(t)" },
        0.0, 50.0, 1000, "ch02_lorenz_time_series.mml");
    
    std::cout << "   Visualized: Lorenz system components as real functions of time\n";
    std::cout << "   See parametric curves demo for the 3D butterfly attractor!\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
/// @brief Entry point for Real Function visualization demos
///////////////////////////////////////////////////////////////////////////////////////////
void Chapter02_RealFunctionVisualization()
{
    std::cout << "\n";
    std::cout << "╔═══════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║           CHAPTER 02: Real Function Visualization                 ║\n";
    std::cout << "╚═══════════════════════════════════════════════════════════════════╝\n";
    
    Chapter02_Demo_ClassicFunctions();
    Chapter02_Demo_FunctionWithDerivatives();
    Chapter02_Demo_LorenzTimeSeries();
}
