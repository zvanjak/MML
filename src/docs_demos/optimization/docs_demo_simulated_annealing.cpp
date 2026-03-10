///////////////////////////////////////////////////////////////////////////////////////////
// docs_demo_simulated_annealing.cpp - Examples of Simulated Annealing with MML
//
// Demonstrates:
// 1. Basic SA optimization on simple functions
// 2. Classic benchmark functions (Rosenbrock, Rastrigin, Ackley, etc.)
// 3. Different cooling schedules (exponential, linear, logarithmic)
// 4. Different neighbor generators (uniform, Gaussian)
// 5. Constrained optimization with bounds
// 6. Higher-dimensional problems
// 7. Comparison with local methods
///////////////////////////////////////////////////////////////////////////////////////////

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"
#include "base/Vector/Vector.h"
#include "SimulatedAnnealing.h"
#include "test_problems/SingleObjective_TestProblems.h"
#endif

#include <iostream>
#include <iomanip>
#include <cmath>
#include <functional>

using namespace MML;
using namespace MML::Optimization;
using namespace MML::TestProblems;

///////////////////////////////////////////////////////////////////////////////////////////
// Helper function to print SA results
///////////////////////////////////////////////////////////////////////////////////////////
void PrintSAResult(const std::string& name, const HeuristicOptimizationResult& result, 
                   const Vector<Real>& trueMin = Vector<Real>(), Real trueMinValue = 0.0)
{
    std::cout << "\n=== " << name << " ===" << std::endl;
    std::cout << std::fixed << std::setprecision(6);
    
    std::cout << "Best solution found: (";
    for (int i = 0; i < result.xbest.size(); ++i) {
        std::cout << result.xbest[i];
        if (i < result.xbest.size() - 1) std::cout << ", ";
    }
    std::cout << ")" << std::endl;
    
    std::cout << "Best objective value: " << result.fbest << std::endl;
    
    if (trueMin.size() > 0) {
        std::cout << "True minimum at: (";
        for (int i = 0; i < trueMin.size(); ++i) {
            std::cout << trueMin[i];
            if (i < trueMin.size() - 1) std::cout << ", ";
        }
        std::cout << ") = " << trueMinValue << std::endl;
        std::cout << "Error in value: " << std::abs(result.fbest - trueMinValue) << std::endl;
    }
    
    std::cout << "Iterations: " << result.iterations << std::endl;
    std::cout << "Function evaluations: " << result.funcEvals << std::endl;
    std::cout << "Accepted moves: " << result.acceptedMoves << std::endl;
    std::cout << "Acceptance rate: " << std::setprecision(2) 
              << (100.0 * result.acceptedMoves / result.funcEvals) << "%" << std::endl;
    std::cout << "Converged: " << (result.converged ? "Yes" : "No") << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
// Example 1: Simple Quadratic Function
///////////////////////////////////////////////////////////////////////////////////////////
void Example_SimpleQuadratic()
{
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "Example 1: Simple Quadratic Function" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    std::cout << "\nMinimize f(x,y) = (x-2)^2 + (y-3)^2" << std::endl;
    std::cout << "Global minimum at (2, 3) with f = 0" << std::endl;
    std::cout << std::endl;
    
    // Define the objective function
    auto quadratic = [](const Vector<Real>& x) -> Real {
        return (x[0] - 2.0) * (x[0] - 2.0) + (x[1] - 3.0) * (x[1] - 3.0);
    };
    
    // Starting point far from minimum
    Vector<Real> x0({10.0, -5.0});
    std::cout << "Starting point: (" << x0[0] << ", " << x0[1] << ")" << std::endl;
    
    // Create SA optimizer with default settings
    SimulatedAnnealing sa(100.0, 1e-10, 5000, 500, 
                          SimulatedAnnealing::StopCriteria::Combined, 42);
    
    // Run optimization
    HeuristicOptimizationResult result = sa.Minimize(quadratic, x0);
    
    Vector<Real> trueMin({2.0, 3.0});
    PrintSAResult("Quadratic Result", result, trueMin, 0.0);
}

///////////////////////////////////////////////////////////////////////////////////////////
// Example 2: Rosenbrock's Banana Function
///////////////////////////////////////////////////////////////////////////////////////////
void Example_Rosenbrock()
{
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "Example 2: Rosenbrock's Banana Function (2D)" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    std::cout << "\nMinimize f(x,y) = (1-x)^2 + 100*(y-x^2)^2" << std::endl;
    std::cout << "Global minimum at (1, 1) with f = 0" << std::endl;
    std::cout << "\nThis is a classic difficult test function with a narrow" << std::endl;
    std::cout << "curved valley. Local methods often struggle here." << std::endl;
    std::cout << std::endl;
    
    // Use centralized Rosenbrock function from TestProblems
    auto rosenbrock = Functions::Rosenbrock();
    
    // Starting point
    Vector<Real> x0({-2.0, 2.0});
    std::cout << "Starting point: (" << x0[0] << ", " << x0[1] << ")" << std::endl;
    
    // SA with more iterations for this difficult problem
    SimulatedAnnealing sa(200.0, 1e-12, 20000, 2000, 
                          SimulatedAnnealing::StopCriteria::Combined, 42);
    
    HeuristicOptimizationResult result = sa.Minimize(rosenbrock, x0);
    
    Vector<Real> trueMin({1.0, 1.0});
    PrintSAResult("Rosenbrock Result", result, trueMin, 0.0);
}

///////////////////////////////////////////////////////////////////////////////////////////
// Example 3: Rastrigin Function (Highly Multimodal)
///////////////////////////////////////////////////////////////////////////////////////////
void Example_Rastrigin()
{
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "Example 3: Rastrigin Function (Highly Multimodal)" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    std::cout << "\nMinimize f(x) = 10n + sum_i [x_i^2 - 10*cos(2*pi*x_i)]" << std::endl;
    std::cout << "Global minimum at origin with f = 0" << std::endl;
    std::cout << "\nThis function has many local minima arranged in a regular" << std::endl;
    std::cout << "pattern - excellent test for global optimization ability." << std::endl;
    std::cout << std::endl;
    
    // Use centralized Rastrigin function from TestProblems
    auto rastrigin = Functions::Rastrigin();
    
    // Start away from origin
    Vector<Real> x0({3.5, -2.5});
    std::cout << "Starting point: (" << x0[0] << ", " << x0[1] << ")" << std::endl;
    std::cout << "f(x0) = " << std::fixed << std::setprecision(4) << rastrigin(x0) << std::endl;
    
    // SA configuration - need good exploration for multimodal
    SimulatedAnnealing sa(500.0, 1e-10, 30000, 3000, 
                          SimulatedAnnealing::StopCriteria::Combined, 42);
    
    HeuristicOptimizationResult result = sa.Minimize(rastrigin, x0);
    
    Vector<Real> trueMin({0.0, 0.0});
    PrintSAResult("Rastrigin Result", result, trueMin, 0.0);
}

///////////////////////////////////////////////////////////////////////////////////////////
// Example 4: Ackley Function
///////////////////////////////////////////////////////////////////////////////////////////
void Example_Ackley()
{
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "Example 4: Ackley Function" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    std::cout << "\nf(x) = -20*exp(-0.2*sqrt(mean(x^2))) - exp(mean(cos(2*pi*x))) + 20 + e" << std::endl;
    std::cout << "Global minimum at origin with f = 0" << std::endl;
    std::cout << "\nAckley has a nearly flat outer region with many local minima" << std::endl;
    std::cout << "and a central global minimum in a deep hole." << std::endl;
    std::cout << std::endl;
    
    // Use centralized Ackley function from TestProblems
    auto ackley = Functions::Ackley();
    
    Vector<Real> x0({4.0, -3.0});
    std::cout << "Starting point: (" << x0[0] << ", " << x0[1] << ")" << std::endl;
    std::cout << "f(x0) = " << std::fixed << std::setprecision(4) << ackley(x0) << std::endl;
    
    SimulatedAnnealing sa(200.0, 1e-10, 15000, 2000, 
                          SimulatedAnnealing::StopCriteria::Combined, 42);
    
    HeuristicOptimizationResult result = sa.Minimize(ackley, x0);
    
    Vector<Real> trueMin({0.0, 0.0});
    PrintSAResult("Ackley Result", result, trueMin, 0.0);
}

///////////////////////////////////////////////////////////////////////////////////////////
// Example 5: Comparing Cooling Schedules
///////////////////////////////////////////////////////////////////////////////////////////
void Example_CoolingSchedules()
{
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "Example 5: Comparing Cooling Schedules" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    std::cout << "\nComparing exponential, linear, and logarithmic cooling" << std::endl;
    std::cout << "on the Sphere function: f(x) = sum(x_i^2)" << std::endl;
    std::cout << std::endl;
    
    // Use centralized Sphere function from TestProblems
    auto sphere = Functions::Sphere();
    
    Vector<Real> x0({5.0, -3.0, 4.0});
    Vector<Real> trueMin({0.0, 0.0, 0.0});
    
    std::cout << "Starting point: (" << x0[0] << ", " << x0[1] << ", " << x0[2] << ")" << std::endl;
    std::cout << "f(x0) = " << sphere(x0) << std::endl;
    std::cout << std::endl;
    
    int maxIter = 10000;
    Real T0 = 100.0;
    Real Tmin = 1e-10;
    
    // 1. Exponential cooling (fast, alpha=0.99)
    std::cout << "--- Exponential Cooling (alpha = 0.99) ---" << std::endl;
    {
        auto cooling = std::make_unique<ExponentialCooling>(0.99);
        auto neighbor = std::make_unique<GaussianNeighborGenerator>(0.5, true, T0, 42);
        
        SimulatedAnnealing sa(std::move(cooling), std::move(neighbor), T0, Tmin, maxIter, 1000,
                              SimulatedAnnealing::StopCriteria::Combined, 42);
        
        HeuristicOptimizationResult result = sa.Minimize(sphere, x0);
        std::cout << "Best value: " << std::fixed << std::setprecision(8) << result.fbest << std::endl;
        std::cout << "Iterations: " << result.iterations << ", Func evals: " << result.funcEvals << std::endl;
    }
    
    // 2. Exponential cooling (slow, alpha=0.999)
    std::cout << "\n--- Exponential Cooling (alpha = 0.999, slower) ---" << std::endl;
    {
        auto cooling = std::make_unique<ExponentialCooling>(0.999);
        auto neighbor = std::make_unique<GaussianNeighborGenerator>(0.5, true, T0, 42);
        
        SimulatedAnnealing sa(std::move(cooling), std::move(neighbor), T0, Tmin, maxIter, 1000,
                              SimulatedAnnealing::StopCriteria::Combined, 42);
        
        HeuristicOptimizationResult result = sa.Minimize(sphere, x0);
        std::cout << "Best value: " << std::fixed << std::setprecision(8) << result.fbest << std::endl;
        std::cout << "Iterations: " << result.iterations << ", Func evals: " << result.funcEvals << std::endl;
    }
    
    // 3. Linear cooling
    std::cout << "\n--- Linear Cooling ---" << std::endl;
    {
        auto cooling = std::make_unique<LinearCooling>(maxIter);
        auto neighbor = std::make_unique<GaussianNeighborGenerator>(0.5, true, T0, 42);
        
        SimulatedAnnealing sa(std::move(cooling), std::move(neighbor), T0, Tmin, maxIter, 1000,
                              SimulatedAnnealing::StopCriteria::Combined, 42);
        
        HeuristicOptimizationResult result = sa.Minimize(sphere, x0);
        std::cout << "Best value: " << std::fixed << std::setprecision(8) << result.fbest << std::endl;
        std::cout << "Iterations: " << result.iterations << ", Func evals: " << result.funcEvals << std::endl;
    }
    
    // 4. Logarithmic cooling (very slow but thorough)
    std::cout << "\n--- Logarithmic Cooling (c = 1.0) ---" << std::endl;
    {
        auto cooling = std::make_unique<LogarithmicCooling>(1.0);
        auto neighbor = std::make_unique<GaussianNeighborGenerator>(0.5, true, T0, 42);
        
        SimulatedAnnealing sa(std::move(cooling), std::move(neighbor), T0, Tmin, maxIter, 1000,
                              SimulatedAnnealing::StopCriteria::Combined, 42);
        
        HeuristicOptimizationResult result = sa.Minimize(sphere, x0);
        std::cout << "Best value: " << std::fixed << std::setprecision(8) << result.fbest << std::endl;
        std::cout << "Iterations: " << result.iterations << ", Func evals: " << result.funcEvals << std::endl;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
// Example 6: Comparing Neighbor Generators
///////////////////////////////////////////////////////////////////////////////////////////
void Example_NeighborGenerators()
{
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "Example 6: Comparing Neighbor Generators" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    std::cout << "\nComparing Uniform vs Gaussian neighbor generation" << std::endl;
    std::cout << "on the Booth function: f(x,y) = (x + 2y - 7)^2 + (2x + y - 5)^2" << std::endl;
    std::cout << "Global minimum at (1, 3) with f = 0" << std::endl;
    std::cout << std::endl;
    
    // Use centralized Booth function from TestProblems
    auto booth = Functions::Booth();
    
    Vector<Real> x0({-5.0, 5.0});
    Vector<Real> trueMin({1.0, 3.0});
    Real T0 = 100.0;
    int maxIter = 8000;
    
    std::cout << "Starting point: (" << x0[0] << ", " << x0[1] << ")" << std::endl;
    std::cout << "f(x0) = " << booth(x0) << std::endl;
    std::cout << std::endl;
    
    // 1. Uniform neighbor generator
    std::cout << "--- Uniform Neighbor Generator (delta = 1.0) ---" << std::endl;
    {
        auto cooling = std::make_unique<ExponentialCooling>(0.995);
        auto neighbor = std::make_unique<UniformNeighborGenerator>(1.0, true, T0, 42);
        
        SimulatedAnnealing sa(std::move(cooling), std::move(neighbor), T0, 1e-10, maxIter, 1000,
                              SimulatedAnnealing::StopCriteria::Combined, 42);
        
        HeuristicOptimizationResult result = sa.Minimize(booth, x0);
        std::cout << "Best: (" << result.xbest[0] << ", " << result.xbest[1] << ")" << std::endl;
        std::cout << "Best value: " << std::fixed << std::setprecision(8) << result.fbest << std::endl;
        std::cout << "Iterations: " << result.iterations << std::endl;
    }
    
    // 2. Gaussian neighbor generator
    std::cout << "\n--- Gaussian Neighbor Generator (sigma = 1.0) ---" << std::endl;
    {
        auto cooling = std::make_unique<ExponentialCooling>(0.995);
        auto neighbor = std::make_unique<GaussianNeighborGenerator>(1.0, true, T0, 42);
        
        SimulatedAnnealing sa(std::move(cooling), std::move(neighbor), T0, 1e-10, maxIter, 1000,
                              SimulatedAnnealing::StopCriteria::Combined, 42);
        
        HeuristicOptimizationResult result = sa.Minimize(booth, x0);
        std::cout << "Best: (" << result.xbest[0] << ", " << result.xbest[1] << ")" << std::endl;
        std::cout << "Best value: " << std::fixed << std::setprecision(8) << result.fbest << std::endl;
        std::cout << "Iterations: " << result.iterations << std::endl;
    }
    
    // 3. Gaussian with temperature scaling
    std::cout << "\n--- Gaussian with Temperature Scaling ---" << std::endl;
    {
        auto cooling = std::make_unique<ExponentialCooling>(0.995);
        auto neighbor = std::make_unique<GaussianNeighborGenerator>(2.0, true, T0, 42);
        
        SimulatedAnnealing sa(std::move(cooling), std::move(neighbor), T0, 1e-10, maxIter, 1000,
                              SimulatedAnnealing::StopCriteria::Combined, 42);
        
        HeuristicOptimizationResult result = sa.Minimize(booth, x0);
        std::cout << "Best: (" << result.xbest[0] << ", " << result.xbest[1] << ")" << std::endl;
        std::cout << "Best value: " << std::fixed << std::setprecision(8) << result.fbest << std::endl;
        std::cout << "Iterations: " << result.iterations << std::endl;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
// Example 7: Bounded Optimization
///////////////////////////////////////////////////////////////////////////////////////////
void Example_BoundedOptimization()
{
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "Example 7: Bounded Optimization" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    std::cout << "\nMinimize Himmelblau's function with bounds" << std::endl;
    std::cout << "f(x,y) = (x^2 + y - 11)^2 + (x + y^2 - 7)^2" << std::endl;
    std::cout << "Has 4 global minima at f = 0" << std::endl;
    std::cout << "Bounds: -5 <= x,y <= 5" << std::endl;
    std::cout << std::endl;
    
    // Use centralized Himmelblau function from TestProblems
    auto himmelblau = Functions::Himmelblau();
    
    Vector<Real> x0({0.0, 0.0});
    Vector<Real> lower({-5.0, -5.0});
    Vector<Real> upper({5.0, 5.0});
    Real T0 = 100.0;
    
    std::cout << "Starting point: (" << x0[0] << ", " << x0[1] << ")" << std::endl;
    std::cout << "f(x0) = " << himmelblau(x0) << std::endl;
    std::cout << std::endl;
    
    // Create bounded neighbor generator
    auto cooling = std::make_unique<ExponentialCooling>(0.995);
    auto neighbor = std::make_unique<GaussianNeighborGenerator>(1.0, true, T0, 42);
    neighbor->SetBounds(lower, upper);
    
    SimulatedAnnealing sa(std::move(cooling), std::move(neighbor), T0, 1e-10, 15000, 2000,
                          SimulatedAnnealing::StopCriteria::Combined, 42);
    
    HeuristicOptimizationResult result = sa.Minimize(himmelblau, x0);
    
    std::cout << "=== Bounded Optimization Result ===" << std::endl;
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Best solution: (" << result.xbest[0] << ", " << result.xbest[1] << ")" << std::endl;
    std::cout << "Best value: " << result.fbest << std::endl;
    std::cout << "Iterations: " << result.iterations << std::endl;
    
    std::cout << "\nKnown global minima (all with f = 0):" << std::endl;
    std::cout << "  (3.0, 2.0)" << std::endl;
    std::cout << "  (-2.805118, 3.131312)" << std::endl;
    std::cout << "  (-3.779310, -3.283186)" << std::endl;
    std::cout << "  (3.584428, -1.848126)" << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
// Example 8: Higher Dimensional Problem (5D Sphere)
///////////////////////////////////////////////////////////////////////////////////////////
void Example_HigherDimensional()
{
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "Example 8: Higher Dimensional Problem (5D)" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    std::cout << "\nMinimize 5D Sphere function: f(x) = sum(x_i^2)" << std::endl;
    std::cout << "Global minimum at origin with f = 0" << std::endl;
    std::cout << std::endl;
    
    // Use centralized Sphere function from TestProblems (works for any dimension)
    auto sphere5d = Functions::Sphere();
    
    Vector<Real> x0({5.0, -3.0, 4.0, -2.0, 6.0});
    Vector<Real> trueMin({0.0, 0.0, 0.0, 0.0, 0.0});
    
    std::cout << "Starting point: (";
    for (int i = 0; i < 5; ++i) {
        std::cout << x0[i];
        if (i < 4) std::cout << ", ";
    }
    std::cout << ")" << std::endl;
    std::cout << "f(x0) = " << sphere5d(x0) << std::endl;
    
    // SA for higher dimensions - need more iterations
    SimulatedAnnealing sa(200.0, 1e-12, 50000, 5000, 
                          SimulatedAnnealing::StopCriteria::Combined, 42);
    
    HeuristicOptimizationResult result = sa.Minimize(sphere5d, x0);
    
    PrintSAResult("5D Sphere Result", result, trueMin, 0.0);
}

///////////////////////////////////////////////////////////////////////////////////////////
// Example 9: Styblinski-Tang Function
///////////////////////////////////////////////////////////////////////////////////////////
void Example_StyblinskiTang()
{
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "Example 9: Styblinski-Tang Function" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    std::cout << "\nf(x) = 0.5 * sum(x_i^4 - 16*x_i^2 + 5*x_i)" << std::endl;
    std::cout << "Global minimum at x_i = -2.903534... with f/n = -39.166166..." << std::endl;
    std::cout << "Testing in 3D (n=3), so f_min = -117.498..." << std::endl;
    std::cout << std::endl;
    
    // Use centralized Styblinski-Tang function from TestProblems
    auto styblinski = Functions::StyblinskiTang();
    
    Vector<Real> x0({2.0, -1.0, 3.0});
    Real trueMinVal = -39.16616570377 * 3;  // Per-dimension minimum * n
    Vector<Real> trueMin({-2.903534, -2.903534, -2.903534});
    
    std::cout << "Starting point: (" << x0[0] << ", " << x0[1] << ", " << x0[2] << ")" << std::endl;
    std::cout << "f(x0) = " << std::fixed << std::setprecision(4) << styblinski(x0) << std::endl;
    
    SimulatedAnnealing sa(200.0, 1e-10, 20000, 2000, 
                          SimulatedAnnealing::StopCriteria::Combined, 42);
    
    HeuristicOptimizationResult result = sa.Minimize(styblinski, x0);
    
    PrintSAResult("Styblinski-Tang Result", result, trueMin, trueMinVal);
}

///////////////////////////////////////////////////////////////////////////////////////////
// Example 10: Schwefel Function (Deceptive)
///////////////////////////////////////////////////////////////////////////////////////////
void Example_Schwefel()
{
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "Example 10: Schwefel Function (Deceptive)" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    std::cout << "\nf(x) = 418.9829*n - sum(x_i * sin(sqrt(|x_i|)))" << std::endl;
    std::cout << "Global minimum at x_i = 420.9687... with f = 0" << std::endl;
    std::cout << "\nThis function is deceptive - the global minimum is" << std::endl;
    std::cout << "geometrically far from the next-best local minima." << std::endl;
    std::cout << "Domain: [-500, 500]^n" << std::endl;
    std::cout << std::endl;
    
    // Use centralized Schwefel function from TestProblems
    auto schwefel = Functions::Schwefel();
    
    Vector<Real> x0({100.0, -200.0});
    Vector<Real> lower({-500.0, -500.0});
    Vector<Real> upper({500.0, 500.0});
    Real T0 = 500.0;
    
    std::cout << "Starting point: (" << x0[0] << ", " << x0[1] << ")" << std::endl;
    std::cout << "f(x0) = " << std::fixed << std::setprecision(4) << schwefel(x0) << std::endl;
    
    // Need high temperature and many iterations for this deceptive function
    auto cooling = std::make_unique<ExponentialCooling>(0.9995);
    auto neighbor = std::make_unique<GaussianNeighborGenerator>(50.0, true, T0, 42);
    neighbor->SetBounds(lower, upper);
    
    SimulatedAnnealing sa(std::move(cooling), std::move(neighbor), T0, 1e-8, 50000, 5000,
                          SimulatedAnnealing::StopCriteria::Combined, 42);
    
    HeuristicOptimizationResult result = sa.Minimize(schwefel, x0);
    
    Vector<Real> trueMin({420.9687, 420.9687});
    PrintSAResult("Schwefel Result", result, trueMin, 0.0);
}

///////////////////////////////////////////////////////////////////////////////////////////
// Example 11: Portfolio Optimization (Financial Application)
///////////////////////////////////////////////////////////////////////////////////////////
void Example_PortfolioOptimization()
{
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "Example 11: Portfolio Optimization" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    std::cout << "\nMinimize portfolio variance given expected return constraint" << std::endl;
    std::cout << "3 assets with returns: 10%, 12%, 14%" << std::endl;
    std::cout << "Covariance matrix and target return = 11%" << std::endl;
    std::cout << std::endl;
    
    // Expected returns
    Vector<Real> returns({0.10, 0.12, 0.14});
    
    // Covariance matrix (simplified)
    // Asset 1: low risk, Asset 2: medium, Asset 3: high
    Real cov[3][3] = {
        {0.04, 0.01, 0.02},
        {0.01, 0.09, 0.03},
        {0.02, 0.03, 0.16}
    };
    
    Real targetReturn = 0.11;
    Real penaltyWeight = 100.0;
    
    // Objective: minimize variance + penalties for constraints
    auto portfolio = [&](const Vector<Real>& w) -> Real {
        // Portfolio variance: w' * Cov * w
        Real variance = 0.0;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                variance += w[i] * cov[i][j] * w[j];
            }
        }
        
        // Penalty for not summing to 1
        Real sumW = w[0] + w[1] + w[2];
        Real penalty1 = penaltyWeight * (sumW - 1.0) * (sumW - 1.0);
        
        // Penalty for not meeting return target
        Real portReturn = w[0] * returns[0] + w[1] * returns[1] + w[2] * returns[2];
        Real penalty2 = penaltyWeight * (portReturn - targetReturn) * (portReturn - targetReturn);
        
        // Penalty for negative weights (no short selling)
        Real penalty3 = 0.0;
        for (int i = 0; i < 3; ++i) {
            if (w[i] < 0.0) penalty3 += penaltyWeight * w[i] * w[i];
        }
        
        return variance + penalty1 + penalty2 + penalty3;
    };
    
    // Start with equal weights
    Vector<Real> w0({0.33, 0.33, 0.34});
    
    std::cout << "Initial equal-weight portfolio:" << std::endl;
    std::cout << "  Weights: (" << w0[0] << ", " << w0[1] << ", " << w0[2] << ")" << std::endl;
    
    // SA with bounds
    Vector<Real> lower({0.0, 0.0, 0.0});
    Vector<Real> upper({1.0, 1.0, 1.0});
    Real T0 = 10.0;
    
    auto cooling = std::make_unique<ExponentialCooling>(0.995);
    auto neighbor = std::make_unique<GaussianNeighborGenerator>(0.05, true, T0, 42);
    neighbor->SetBounds(lower, upper);
    
    SimulatedAnnealing sa(std::move(cooling), std::move(neighbor), T0, 1e-10, 20000, 2000,
                          SimulatedAnnealing::StopCriteria::Combined, 42);
    
    HeuristicOptimizationResult result = sa.Minimize(portfolio, w0);
    
    std::cout << "\n=== Optimized Portfolio ===" << std::endl;
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "Weights: Asset1=" << result.xbest[0] 
              << ", Asset2=" << result.xbest[1] 
              << ", Asset3=" << result.xbest[2] << std::endl;
    
    Real sumW = result.xbest[0] + result.xbest[1] + result.xbest[2];
    std::cout << "Sum of weights: " << sumW << " (should be 1.0)" << std::endl;
    
    Real portReturn = result.xbest[0] * returns[0] + result.xbest[1] * returns[1] + result.xbest[2] * returns[2];
    std::cout << "Expected return: " << (portReturn * 100) << "% (target: " << (targetReturn * 100) << "%)" << std::endl;
    
    // Calculate portfolio variance
    Real variance = 0.0;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            variance += result.xbest[i] * cov[i][j] * result.xbest[j];
        }
    }
    std::cout << "Portfolio variance: " << variance << std::endl;
    std::cout << "Portfolio std dev: " << std::sqrt(variance) << " (" << (std::sqrt(variance) * 100) << "%)" << std::endl;
    std::cout << "Iterations: " << result.iterations << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
// Example 12: Engineering Design - Spring Optimization
///////////////////////////////////////////////////////////////////////////////////////////
void Example_SpringDesign()
{
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "Example 12: Engineering Design - Spring Optimization" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    std::cout << "\nMinimize spring weight: f(d,D,N) = (N+2)*D*d^2" << std::endl;
    std::cout << "  d = wire diameter" << std::endl;
    std::cout << "  D = mean coil diameter" << std::endl;
    std::cout << "  N = number of active coils" << std::endl;
    std::cout << "\nSubject to constraints (handled via penalty):" << std::endl;
    std::cout << "  Deflection, shear stress, surge frequency, outer diameter" << std::endl;
    std::cout << std::endl;
    
    // Simplified spring weight minimization
    // Real problem has complex constraints, we simplify here
    auto springWeight = [](const Vector<Real>& x) -> Real {
        Real d = x[0];  // wire diameter
        Real D = x[1];  // mean coil diameter
        Real N = x[2];  // number of coils
        
        // Objective: minimize weight
        Real weight = (N + 2.0) * D * d * d;
        
        // Simplified constraints as penalties
        Real penalty = 0.0;
        
        // Constraint 1: D/d ratio (spring index) should be 4-12
        Real springIndex = D / d;
        if (springIndex < 4.0) penalty += 1000.0 * (4.0 - springIndex) * (4.0 - springIndex);
        if (springIndex > 12.0) penalty += 1000.0 * (springIndex - 12.0) * (springIndex - 12.0);
        
        // Constraint 2: minimum number of coils
        if (N < 3.0) penalty += 1000.0 * (3.0 - N) * (3.0 - N);
        
        // Constraint 3: outer diameter limit (D + d <= 3)
        if (D + d > 3.0) penalty += 1000.0 * (D + d - 3.0) * (D + d - 3.0);
        
        return weight + penalty;
    };
    
    Vector<Real> x0({0.3, 1.5, 10.0});  // Initial guess
    Vector<Real> lower({0.05, 0.25, 3.0});
    Vector<Real> upper({1.0, 3.0, 50.0});
    Real T0 = 50.0;
    
    std::cout << "Initial design: d=" << x0[0] << ", D=" << x0[1] << ", N=" << x0[2] << std::endl;
    
    auto cooling = std::make_unique<ExponentialCooling>(0.995);
    auto neighbor = std::make_unique<GaussianNeighborGenerator>(0.2, true, T0, 42);
    neighbor->SetBounds(lower, upper);
    
    SimulatedAnnealing sa(std::move(cooling), std::move(neighbor), T0, 1e-10, 30000, 3000,
                          SimulatedAnnealing::StopCriteria::Combined, 42);
    
    HeuristicOptimizationResult result = sa.Minimize(springWeight, x0);
    
    std::cout << "\n=== Optimized Spring Design ===" << std::endl;
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "Wire diameter (d): " << result.xbest[0] << std::endl;
    std::cout << "Mean coil diameter (D): " << result.xbest[1] << std::endl;
    std::cout << "Number of coils (N): " << result.xbest[2] << std::endl;
    std::cout << "Spring index (D/d): " << result.xbest[1] / result.xbest[0] << std::endl;
    std::cout << "Outer diameter (D+d): " << result.xbest[1] + result.xbest[0] << std::endl;
    
    Real finalWeight = (result.xbest[2] + 2.0) * result.xbest[1] * result.xbest[0] * result.xbest[0];
    std::cout << "Spring weight: " << finalWeight << std::endl;
    std::cout << "Iterations: " << result.iterations << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
// Main function - run all examples
///////////////////////////////////////////////////////////////////////////////////////////
void Docs_Demo_SimulatedAnnealing()
{
    std::cout << std::endl;
    std::cout << "###########################################################################" << std::endl;
    std::cout << "#                                                                         #" << std::endl;
    std::cout << "#                   SIMULATED ANNEALING EXAMPLES                          #" << std::endl;
    std::cout << "#                                                                         #" << std::endl;
    std::cout << "###########################################################################" << std::endl;
    
    // Basic examples
    Example_SimpleQuadratic();
    Example_Rosenbrock();
    
    // Multimodal benchmark functions
    Example_Rastrigin();
    Example_Ackley();
    
    // Algorithm configuration
    Example_CoolingSchedules();
    Example_NeighborGenerators();
    
    // Constrained and bounded
    Example_BoundedOptimization();
    
    // Higher dimensions
    Example_HigherDimensional();
    
    // More benchmarks
    Example_StyblinskiTang();
    Example_Schwefel();
    
    // Real-world applications
    Example_PortfolioOptimization();
    Example_SpringDesign();
    
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "All Simulated Annealing examples completed!" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
}
