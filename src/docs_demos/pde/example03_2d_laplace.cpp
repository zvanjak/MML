/**
 * Example 3: 2D Laplace Equation
 * 
 * Solves: ∇²u = 0 on [0,1] × [0,1]
 * with boundary conditions:
 *   u(x,0) = 0           (bottom)
 *   u(x,1) = sin(πx)     (top)
 *   u(0,y) = 0           (left)
 *   u(1,y) = 0           (right)
 * 
 * This represents a steady-state temperature distribution
 * in a square plate with heated top edge.
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
#include "tools/Visualizer.h"

using namespace MML;
using namespace MML::PDE;

void Docs_Demo_PDE_2D_Laplace() {
    // Spatial discretization - use finer grid for better visualization
    MML::PDE::Rectangle<double> domain(0.0, 1.0, 0.0, 1.0);
    Grid2D<double> grid(domain, 80, 80);  // 80x80 cells = 81x81 nodes
    
    std::cout << "2D Laplace Equation Results:\n";
    std::cout << "============================\n\n";
    
    std::cout << "Parameters:\n";
    std::cout << "  Domain: [0,1] × [0,1]\n";
    std::cout << "  Grid: " << grid.numNodesX() << " × " << grid.numNodesY() << " nodes\n\n";
    
    // Set boundary conditions
    BoundaryConditions2D<double> bc;
    bc.setBottom(BoundaryCondition2D<double>::Dirichlet(0.0));
    bc.setTop(BoundaryCondition2D<double>::Dirichlet(
        std::function<double(double, double)>([](double x, double y) { 
            return std::sin(M_PI * x); 
        })
    ));
    bc.setLeft(BoundaryCondition2D<double>::Dirichlet(0.0));
    bc.setRight(BoundaryCondition2D<double>::Dirichlet(0.0));
    
    // Create Poisson solver (Laplace is Poisson with f=0)
    PoissonSolver2D<double> laplace(grid, bc);
    laplace.setSource([](double x, double y) { return 0.0; });  // Laplace: ∇²u = 0
    
    // Solve
    std::cout << "Solving Laplace equation...\n";
    auto solution = laplace.solve();
    std::cout << "Solution computed!\n\n";
    
    // Print results at selected interior points
    std::cout << "Solution at selected points:\n";
    std::cout << "(x, y)\t\t\tTemperature\n";
    std::cout << "----------------------------------------\n";
    
    int nx = grid.numNodesX();
    int ny = grid.numNodesY();
    
    const std::vector<std::pair<int, int>> sample_points = {
        {nx/4, ny/4},    // Near bottom-left
        {nx/2, ny/4},    // Bottom center
        {3*nx/4, ny/4},  // Near bottom-right
        {nx/2, ny/2},    // Center
        {nx/2, 3*ny/4}   // Top center
    };
    
    for (const auto& [i, j] : sample_points) {
        double x = grid.x(i);
        double y = grid.y(j);
        double value = solution(i, j);
        
        std::cout << "(" << x << ", " << y << ")\t\t" << value << "\n";
    }
    
    // Verify boundary conditions
    std::cout << "\nBoundary verification:\n";
    std::cout << "  Top center: " << solution(nx/2, ny-1) 
              << " (expected: " << std::sin(M_PI * 0.5) << ")\n";
    std::cout << "  Bottom center: " << solution(nx/2, 0) << " (expected: 0)\n";
    std::cout << "  Left center: " << solution(0, ny/2) << " (expected: 0)\n";
    std::cout << "  Right center: " << solution(nx-1, ny/2) << " (expected: 0)\n";
    
    // ====================================================================
    // VISUALIZATION: First PDE solution visualization! 🎉
    // Scale domain by 100 (0-1 → 0-100) and values by 100 for better visibility
    // ====================================================================
    std::cout << "\n========================================\n";
    std::cout << "Launching temperature visualization...\n";
    std::cout << "(Domain scaled 100x, values scaled 100x)\n";
    std::cout << "========================================\n";
    
    auto vizResult = Visualizer::VisualizeGridFunction2D(
        solution, 
        "Temperature Distribution - 2D Laplace Equation (Heated Top Edge)",
        "pde_laplace_2d_temperature.mml",
        10.0,   // Scale x,y coordinates by 10
        10.0    // Scale temperature values by 10
    );
    
    if (!vizResult.success) {
        std::cerr << "Visualization failed: " << vizResult.errorMessage << "\n";
        return;
    }
    
    std::cout << "Visualization complete!\n";
}
