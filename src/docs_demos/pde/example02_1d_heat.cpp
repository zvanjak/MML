/**
 * Example 2: 1D Heat Equation
 * 
 * Solves: ∂u/∂t = α ∂²u/∂x² on (0,1) × (0,T]
 * with initial condition u(x,0) = sin(πx)
 * and boundary conditions u(0,t) = u(1,t) = 0
 * 
 * Analytical solution: u(x,t) = sin(πx)·exp(-π²αt)
 * 
 * From: docs/pde/Quick_Start_Guide.md
 */

#include <iostream>
#include <cmath>

// Define M_PI for MSVC
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "MMLBase.h"
#include "parabolic/HeatSolver.h"

using namespace MML::PDE;

// Initial condition: sin(πx)
namespace {
double ex2_initial_condition(double x) {
    return std::sin(M_PI * x);
}

// Analytical solution for comparison
double ex2_analytical_solution(double x, double t, double alpha) {
    return std::sin(M_PI * x) * std::exp(-M_PI * M_PI * alpha * t);
}
}

void Docs_Demo_PDE_1D_Heat() {
    // Physical parameters
    const double alpha = 0.01;      // Thermal diffusivity
    const double T_final = 0.5;     // Final time
    const double dt = 0.0001;       // Time step
    
    // Spatial discretization
    Interval<double> domain(0.0, 1.0);
    Grid1D<double> grid(domain, 100);  // 100 cells = 101 nodes
    
    std::cout << "1D Heat Equation Results:\n";
    std::cout << "=========================\n\n";
    
    std::cout << "Parameters:\n";
    std::cout << "  α = " << alpha << "\n";
    std::cout << "  T_final = " << T_final << "\n";
    std::cout << "  dt = " << dt << "\n";
    std::cout << "  Grid: " << grid.numNodes() << " nodes\n\n";
    
    // Set boundary conditions (homogeneous Dirichlet)
    auto bc = homogeneousDirichlet1D<double>();
    
    // Create solver
    HeatSolver1D<double> heat_solver(grid, bc, alpha);
    heat_solver.setInitialCondition(ex2_initial_condition);
    
    // Solve using Crank-Nicolson (unconditionally stable, 2nd order accurate)
    auto result = heat_solver.solve(T_final, dt, TimeScheme::CrankNicolson);
    
    // Get final solution from solver
    const auto& solution = heat_solver.getSolution();
    
    // Print results
    std::cout << "Solver completed " << result.steps << " time steps\n";
    std::cout << "Solution stable: " << (result.stable ? "Yes" : "No") << "\n\n";
    
    std::cout << "Solution at t = " << T_final << ":\n";
    std::cout << "x\t\tNumerical\tAnalytical\tError\n";
    std::cout << "--------------------------------------------------------\n";
    
    const std::vector<int> sample_indices = {0, 25, 50, 75, 100};
    double max_error = 0.0;
    
    for (int i : sample_indices) {
        double x = grid.x(i);
        double numerical = solution[i];
        double analytical = ex2_analytical_solution(x, T_final, alpha);
        double error = std::abs(numerical - analytical);
        max_error = std::max(max_error, error);
        
        std::cout << x << "\t\t" 
                  << numerical << "\t" 
                  << analytical << "\t"
                  << error << "\n";
    }
    
    std::cout << "\nMaximum Error: " << max_error << "\n";
    std::cout << "\nExpected decay: exp(-π²αT) = " 
              << std::exp(-M_PI * M_PI * alpha * T_final) << "\n";
}
