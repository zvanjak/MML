///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        readme10_root_finding.cpp                                           ///
///  Description: README example - Root Finding Algorithms                            ///
///               Demonstrates Bisection, Newton, Brent, Ridders methods              ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///////////////////////////////////////////////////////////////////////////////////////////

#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/MMLBase.h>

#include <mml/base/Function.h>
#include <mml/algorithms/RootFinding.h>
#endif

#include <iostream>
#include <iomanip>

using namespace MML;

void Readme_RootFinding()
{
    std::cout << std::endl;
    std::cout << "=== Root Finding Algorithms ===" << std::endl;

    // Find root of f(x) = x³ - 2x - 5 (has root near x ≈ 2.0945)
    RealFunction f{[](Real x) { return x*x*x - 2*x - 5; }};

    std::cout << "Finding root of f(x) = x³ - 2x - 5 in [2, 3]" << std::endl;
    std::cout << std::setprecision(15);

    // Compare different methods
    Real root_bisect = RootFinding::FindRootBisection(f, 2.0, 3.0, 1e-12);
    Real root_brent  = RootFinding::FindRootBrent(f, 2.0, 3.0, 1e-12);
    Real root_newton = RootFinding::FindRootNewton(f, 2.0, 3.0, 1e-12);
    Real root_ridder = RootFinding::FindRootRidders(f, 2.0, 3.0, 1e-12);

    std::cout << std::endl << "Method comparison:" << std::endl;
    std::cout << "  Bisection: " << root_bisect << std::endl;
    std::cout << "  Brent:     " << root_brent << std::endl;
    std::cout << "  Newton:    " << root_newton << std::endl;
    std::cout << "  Ridders:   " << root_ridder << std::endl;

    // Get detailed convergence info using config
    RootFinding::RootFindingConfig config;
    config.x_tolerance = 1e-14;
    config.f_tolerance = 1e-14;
    config.max_iterations = 100;
    
    auto result = RootFinding::FindRootBrent(f, 2.0, 3.0, config);
    std::cout << std::endl << "Brent with detailed config:" << std::endl;
    std::cout << "  Root:       " << result.root << std::endl;
    std::cout << "  Iterations: " << result.iterations_used << std::endl;
    std::cout << "  f(root) =   " << result.function_value << std::endl;
    std::cout << "  Converged:  " << (result.converged ? "yes" : "no") << std::endl;

    // Find and refine every detectable real root, including tangent roots
    auto allRoots = RootFinding::FindAllRealRootsInInterval(
        [](Real x) { return (x - 1.0) * (x - 2.0) * (x - 2.0) * (x - 3.0); },
        0.0, 4.0);
    
    std::cout << std::endl << "All roots with multiplicity collapsed:" << std::endl;
    for (Real root : allRoots.roots)
        std::cout << "    x = " << root << std::endl;

    // Solve a nonlinear system F(x)=0 with the shared numerical Jacobian
    auto system = [](const Vector<Real>& point) {
        return Vector<Real>({point[0] * point[0] + point[1] * point[1] - Real{1},
                             point[0] - point[1]});
    };
    auto systemResult = RootFinding::SolveNonlinearSystemNewton(
        system, Vector<Real>({0.8, 0.4}));
    std::cout << "Nonlinear solution: " << systemResult.solution[0]
              << ", " << systemResult.solution[1] << std::endl;

    std::cout << std::endl;
}
