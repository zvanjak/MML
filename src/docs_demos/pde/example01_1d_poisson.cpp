/**
 * Example 1: 1D Poisson Equation
 * 
 * Solves: -u'' = sin(πx) on (0,1)
 * with u(0) = 0, u(1) = 0
 * 
 * Analytical solution: u(x) = sin(πx)/π²
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
#include "elliptic/PoissonSolver.h"

using namespace MML::PDE;

// Right-hand side function: sin(πx)
namespace {
double source_function(double x) {
    return std::sin(M_PI * x);
}
}

// Analytical solution for comparison
namespace {
double analytical_solution(double x) {
    return std::sin(M_PI * x) / (M_PI * M_PI);
}
}

void Docs_Demo_PDE_1D_Poisson() {
    // Step 1: Create grid
    Interval<double> domain(0.0, 1.0);
    Grid1D<double> grid(domain, 100);  // 100 cells = 101 nodes
    
    // Step 2: Set boundary conditions (homogeneous Dirichlet)
    auto bc = homogeneousDirichlet1D<double>();
    
    // Step 3: Create Poisson solver
    PoissonSolver1D<double> poisson(grid, bc);
    poisson.setSource(source_function);
    
    // Step 4: Solve
    auto solution = poisson.solve();
    
    // Print results
    std::cout << "1D Poisson Equation Results:\n";
    std::cout << "============================\n\n";
    
    std::cout << "Selected points comparison:\n";
    std::cout << "x\t\tNumerical\tAnalytical\tError\n";
    std::cout << "--------------------------------------------------------\n";
    
    const std::vector<int> sample_indices = {0, 25, 50, 75, 100};
    double max_error = 0.0;
    
    for (int i : sample_indices) {
        double x = grid.x(i);
        double numerical = solution(i);
        double analytical = analytical_solution(x);
        double error = std::abs(numerical - analytical);
        max_error = std::max(max_error, error);
        
        std::cout << x << "\t\t" 
                  << numerical << "\t" 
                  << analytical << "\t"
                  << error << "\n";
    }
    
    std::cout << "\nMaximum Error: " << max_error << "\n";
}
