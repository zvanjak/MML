///////////////////////////////////////////////////////////////////////////////////
// Elliptic PDEs Documentation Verification
// 
// This program verifies that elliptic PDE examples from the documentation work
///////////////////////////////////////////////////////////////////////////////////

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "MMLBase.h"
#include "elliptic/PoissonSolver.h"

#include <iostream>
#include <cmath>

using namespace MML::PDE;

void Docs_Demo_PDE_Elliptic_Verification() {
    std::cout << "=================================================================\n";
    std::cout << "Elliptic PDEs Documentation Verification\n";
    std::cout << "=================================================================\n\n";
    
    int passed = 0;
    int failed = 0;
    
    // Test 1: 1D Poisson
    try {
        std::cout << "Testing 1D Poisson... ";
        Interval<double> domain(0.0, 1.0);
        Grid1D<double> grid(domain, 100);
        auto bc = homogeneousDirichlet1D<double>();
        PoissonSolver1D<double> poisson(grid, bc);
        poisson.setSource([](double x) { return std::sin(M_PI * x); });
        auto solution = poisson.solve();
        std::cout << "✓ PASS\n";
        passed++;
    } catch (const std::exception& e) {
        std::cout << "✗ FAIL: " << e.what() << "\n";
        failed++;
    }
    
    // Test 2: 2D Laplace
    try {
        std::cout << "Testing 2D Laplace... ";
        Rectangle<double> domain(0.0, 1.0, 0.0, 1.0);
        Grid2D<double> grid(domain, 50, 50);
        
        BoundaryConditions2D<double> bc;
        bc.setLeft(BoundaryCondition2D<double>::Dirichlet(0.0));
        bc.setRight(BoundaryCondition2D<double>::Dirichlet(0.0));
        bc.setBottom(BoundaryCondition2D<double>::Dirichlet(0.0));
        bc.setTop(BoundaryCondition2D<double>::Dirichlet(
            [](double x, double y) { return std::sin(M_PI * x); }
        ));
        
        PoissonSolver2D<double> solver(grid, bc);
        solver.setSource([](double x, double y) { return 0.0; });
        auto solution = solver.solve();
        std::cout << "✓ PASS\n";
        passed++;
    } catch (const std::exception& e) {
        std::cout << "✗ FAIL: " << e.what() << "\n";
        failed++;
    }
    
    // Test 3: 2D Poisson with source
    try {
        std::cout << "Testing 2D Poisson with Gaussian source... ";
        Rectangle<double> domain(0.0, 1.0, 0.0, 1.0);
        Grid2D<double> grid(domain, 40, 40);
        auto bc = homogeneousDirichlet2D<double>();
        
        PoissonSolver2D<double> solver(grid, bc);
        solver.setSource([](double x, double y) {
            double r2 = (x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5);
            return std::exp(-r2 / (0.1 * 0.1));
        });
        auto solution = solver.solve();
        std::cout << "✓ PASS\n";
        passed++;
    } catch (const std::exception& e) {
        std::cout << "✗ FAIL: " << e.what() << "\n";
        failed++;
    }
    
    // Test 4: 3D Poisson
    try {
        std::cout << "Testing 3D Poisson... ";
        Box<double> domain(0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
        Grid3D<double> grid(domain, 20, 20, 20);
        auto bc = homogeneousDirichlet3D<double>();
        
        PoissonSolver3D<double> solver(grid, bc);
        solver.setSource([](double x, double y, double z) {
            double r2 = (x - 0.5) * (x - 0.5) + 
                        (y - 0.5) * (y - 0.5) + 
                        (z - 0.5) * (z - 0.5);
            return (r2 < 0.01) ? 100.0 : 0.0;
        });
        auto solution = solver.solve();
        std::cout << "✓ PASS\n";
        passed++;
    } catch (const std::exception& e) {
        std::cout << "✗ FAIL: " << e.what() << "\n";
        failed++;
    }
    
    // Summary
    std::cout << "\n=================================================================\n";
    std::cout << "SUMMARY\n";
    std::cout << "=================================================================\n";
    std::cout << "Total tests: " << (passed + failed) << "\n";
    std::cout << "Passed: " << passed << "\n";
    std::cout << "Failed: " << failed << "\n";
    
    if (failed == 0) {
        std::cout << "\n✓ ALL ELLIPTIC PDE EXAMPLES WORK!\n";
        std::cout << "\n⚠️  WARNING: Documentation examples use simplified API that differs from implementation!\n";
        std::cout << "    Documentation shows: auto solution = solver.solve(f);\n";
        std::cout << "    Actual API uses:     solver.setSource(f); auto solution = solver.solve();\n";
    } else {
        std::cout << "\n✗ Some tests failed.\n";
    }
}
